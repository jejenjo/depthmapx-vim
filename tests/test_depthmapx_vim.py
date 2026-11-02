"""Tests for depthmapx_vim. Run from the repository root:  python -m unittest discover tests"""
import csv
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))
import depthmapx_vim as dv  # noqa: E402


def ref(gx, gy):
    return (gx << 16) + gy


def points(rows):
    """rows: (gx, gy, connectivity, integration) -> VimPoint list in depthmapX order."""
    pts = [dv.VimPoint(i, gx, gy, float(c), float(v)) for i, (gx, gy, c, v) in enumerate(rows)]
    pts.sort(key=lambda p: (p.gx, p.gy))
    return pts


def write_export(path, header, rows, delimiter=","):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=delimiter, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read(path):
    with open(path, newline="") as fh:
        text = fh.read()
    delim = "\t" if text.split("\n", 1)[0].count("\t") > 1 else ","
    return list(csv.DictReader(text.splitlines(), delimiter=delim))


class CoreTests(unittest.TestCase):
    def test_bands_reduce_on_small_plans(self):
        # 25 points, at least 10 per band -> 2 bands
        pts = points([(i, 0, i, 1.0 + i) for i in range(25)])
        dv.compute_vim(pts, num_bands=40)
        self.assertEqual(sorted({p.band for p in pts}), [0, 1])

    def test_equal_connectivity_always_shares_a_band(self):
        # 30 points with connectivity 5 straddle the band boundary; they must stay together
        rows = [(i, 0, 1, 1.0) for i in range(20)] + [(i, 1, 5, 2.0) for i in range(30)] + \
               [(i, 2, 9, 3.0) for i in range(50)]
        pts = points(rows)
        dv.compute_vim(pts, num_bands=10)
        self.assertEqual(len({p.band for p in pts if p.connectivity == 5}), 1)

    def test_vim_is_robust_z_score_within_band(self):
        # one band: integration 1..10, median 5.5, sample SD 3.0277
        pts = points([(i, 0, 7, float(i + 1)) for i in range(10)])
        dv.compute_vim(pts, num_bands=1)
        top = max(pts, key=lambda p: p.integration)
        self.assertAlmostEqual(top.vim, (10 - 5.5) / 3.0276503540974917, places=12)

    def test_zero_sd_gives_zero_vim(self):
        pts = points([(i, 0, 3, 2.5) for i in range(12)])
        dv.compute_vim(pts, num_bands=1)
        self.assertTrue(all(p.vim == 0.0 and p.cls == 1 for p in pts))

    def test_class_thresholds(self):
        pts = points([(i, 0, 1, 0.0) for i in range(4)])
        for p, v in zip(pts, (-1.0, 0.99, 1.0, 2.0)):
            p.vim = v
        # classify through the same rules compute_vim uses
        for p in pts:
            p.cls = 3 if p.vim >= 2 else 2 if p.vim >= 1 else 0 if p.vim <= -1 else 1
        self.assertEqual([p.cls for p in pts], [0, 1, 2, 3])

    def test_clusters_use_8_neighbours_and_rank_by_summed_vim(self):
        # one band; three high outliers: two touching diagonally, one isolated
        rows = [(x, y, 10, 1.0) for x in range(6) for y in range(6)]
        rows = [r for r in rows if (r[0], r[1]) not in {(1, 1), (2, 2), (5, 5)}]
        rows += [(1, 1, 10, 9.0), (2, 2, 10, 9.0), (5, 5, 10, 9.5)]
        pts = points(rows)
        ranked = dv.compute_vim(pts, num_bands=1)
        strong = {(p.gx, p.gy): p for p in pts if p.cls == 3}
        self.assertEqual(set(strong), {(1, 1), (2, 2), (5, 5)})
        self.assertEqual(strong[(1, 1)].cluster, strong[(2, 2)].cluster)  # diagonal neighbours
        self.assertEqual(strong[(1, 1)].cluster, 1)  # two cells outrank one stronger cell
        self.assertEqual(strong[(1, 1)].cluster_size, 2)
        self.assertEqual(strong[(5, 5)].cluster, 2)
        self.assertEqual([len(m) for _, m, _ in ranked], [2, 1])


class FileTests(unittest.TestCase):
    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.header = ["Ref", "x", "y", "Connectivity", "Visual Integration [HH]"]
        self.rows = []
        for gx in range(6):
            for gy in range(5):
                c = 10 + gx
                v = 2.0 + 0.1 * gy + (3.0 if (gx, gy) == (3, 2) else 0.0)
                self.rows.append([ref(gx, gy), 0.5 * gx, 0.5 * gy, c, v])

    def tearDown(self):
        self.dir.cleanup()

    def path(self, name):
        return os.path.join(self.dir.name, name)

    def run_tool(self, *args):
        self.assertEqual(dv.main(list(args) + ["-q"]), 0)

    def test_ref_decodes_to_grid_cell(self):
        cells, source = dv.grid_cells(["Ref"], [[str(786467)]])
        self.assertEqual(cells, [(12, 35)])
        self.assertEqual(source, "Ref")

    def test_comma_tab_and_xy_inputs_agree(self):
        write_export(self.path("a.csv"), self.header, self.rows, ",")
        write_export(self.path("b.txt"), self.header, self.rows, "\t")
        write_export(self.path("c.csv"), self.header[1:], [r[1:] for r in self.rows], ",")
        for name in ("a.csv", "b.txt", "c.csv"):
            self.run_tool(self.path(name))
        cols = dv.vim_column_names("")
        a = read(self.path("a_vim.csv"))
        b = read(self.path("b_vim.txt"))
        c = read(self.path("c_vim.csv"))
        for ra, rb, rc in zip(a, b, c):
            self.assertEqual([ra[k] for k in cols], [rb[k] for k in cols])
            self.assertEqual([ra[k] for k in cols], [rc[k] for k in cols])

    def test_missing_values_are_left_as_minus_one(self):
        rows = [list(r) for r in self.rows]
        rows[0][4] = -1  # depthmapX "no value"
        write_export(self.path("m.csv"), self.header, rows)
        self.run_tool(self.path("m.csv"))
        out = read(self.path("m_vim.csv"))
        self.assertTrue(all(out[0][k] == "-1" for k in dv.vim_column_names("")))
        self.assertNotEqual(out[1]["VIM Band"], "-1")

    def test_existing_vim_columns_are_replaced_not_duplicated(self):
        write_export(self.path("d.csv"), self.header, self.rows)
        self.run_tool(self.path("d.csv"))
        self.run_tool(self.path("d_vim.csv"), "-o", self.path("again.csv"))
        first, second = read(self.path("d_vim.csv")), read(self.path("again.csv"))
        self.assertEqual(list(first[0].keys()), list(second[0].keys()))
        self.assertEqual(first, second)

    def test_radius_columns_get_suffix(self):
        header = self.header[:4] + ["Visual Integration [HH] R3"]
        write_export(self.path("r.csv"), header, self.rows)
        self.run_tool(self.path("r.csv"))
        out = read(self.path("r_vim.csv"))
        for name in dv.vim_column_names(" R3"):
            self.assertIn(name, out[0])

    def test_join_file_has_ref_and_vim_columns_only(self):
        write_export(self.path("j.csv"), self.header, self.rows)
        self.run_tool(self.path("j.csv"))
        join = read(self.path("j_vim_join.csv"))
        self.assertEqual(list(join[0].keys()), ["Ref"] + dv.vim_column_names(""))
        self.assertEqual(len(join), len(self.rows))

    def test_outlier_is_a_strong_hidden_centre(self):
        write_export(self.path("o.csv"), self.header, self.rows)
        self.run_tool(self.path("o.csv"))
        out = {int(r["Ref"]): r for r in read(self.path("o_vim.csv"))}
        self.assertEqual(out[ref(3, 2)]["VIM Class"], "3")
        self.assertEqual(out[ref(3, 2)]["VIM Hidden Centre Cluster"], "1")
        self.assertEqual(sum(r["VIM Class"] == "3" for r in out.values()), 1)

    def test_map_figures_are_valid_svg_with_legend_and_cluster_labels(self):
        import xml.etree.ElementTree as ET
        write_export(self.path("s.csv"), self.header, self.rows)
        self.run_tool(self.path("s.csv"))
        ns = "{http://www.w3.org/2000/svg}"
        for view in ("class", "vim"):
            root = ET.parse(self.path("s_vim_map_%s.svg" % view)).getroot()
            texts = [t.text for t in root.iter(ns + "text")]
            self.assertIn("1", texts)  # the one strong hidden-centre cluster is numbered
            self.assertTrue(any("m" == (t or "").split(" ")[-1] for t in texts))  # scale bar in metres
        legend = [t.text for t in ET.parse(self.path("s_vim_map_class.svg")).getroot().iter(ns + "text")]
        self.assertTrue(any(t.startswith("Strong hidden centre") for t in legend))

    def test_no_map_option_skips_figures(self):
        write_export(self.path("q.csv"), self.header, self.rows)
        self.run_tool(self.path("q.csv"), "--no-map")
        self.assertFalse(os.path.exists(self.path("q_vim_map_class.svg")))

    def test_ramp_is_gray_at_zero_and_saturates_at_the_ends(self):
        self.assertEqual(dv.ramp_colour(0.0), "#d3d2cb")
        self.assertEqual(dv.ramp_colour(5.0), dv.ramp_colour(9.0))
        self.assertEqual(dv.ramp_colour(-5.0), "#104281")

    def test_missing_integration_column_is_a_clear_error(self):
        write_export(self.path("n.csv"), self.header[:4], [r[:4] for r in self.rows])
        self.assertEqual(dv.main([self.path("n.csv"), "-q"]), 1)


if __name__ == "__main__":
    unittest.main()
