# SPDX-FileCopyrightText: Contributors to PyPSA-DE
# SPDX-License-Identifier: MIT

from pathlib import Path

import numpy as np
import pypsa
import pytest

from scripts.add_transmission_projects_and_dlr import attach_transmission_projects
from scripts.build_transmission_projects import add_line_types

LINE_TYPES = Path(__file__).parents[1] / "data/transmission_projects/line_types.csv"


def test_reviewed_line_types_extend_catalogue_once():
    n = pypsa.Network()
    added = add_line_types(n, LINE_TYPES)
    assert n.line_types.loc[added.index, "i_nom"].tolist() == [3.6, 2.36]
    with pytest.raises(ValueError, match="already exist"):
        add_line_types(n, LINE_TYPES)


def test_exported_line_types_are_registered_before_lines(tmp_path):
    add_line_types(pypsa.Network(), LINE_TYPES).to_csv(tmp_path / "line_types.csv")
    (tmp_path / "new_lines.csv").write_text(
        ",bus0,bus1,type,length,num_parallel,s_nom\n"
        "l,a,b,HTLS 2-bundle 220.0,2,2,1798.7\n"
    )
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], v_nom=220.0)
    attach_transmission_projects(
        n, [tmp_path / "line_types.csv", tmp_path / "new_lines.csv"]
    )
    n.calculate_dependent_values()
    assert np.isclose(n.lines.at["l", "x"], 0.302 * 2 / 2)


def test_split_pieces_replace_their_base_branch(tmp_path):
    (tmp_path / "new_buses.csv").write_text(",v_nom\nj,380\n")
    (tmp_path / "new_lines.csv").write_text(
        ",bus0,bus1,length,split_of\nbase#0,a,j,1,base\nbase#1,j,b,1,base\nnew,a,b,2,\n"
    )
    n = pypsa.Network()
    n.add("Bus", ["a", "b"], v_nom=380.0)
    n.add("Line", "base", bus0="a", bus1="b", length=2)
    attach_transmission_projects(
        n, [tmp_path / "new_buses.csv", tmp_path / "new_lines.csv"]
    )
    assert sorted(n.lines.index) == ["base#0", "base#1", "new"]
