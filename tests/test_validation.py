import tempfile
import unittest
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from src.validation import (
    validate_embedding_table,
    validate_full_graph_sidecar,
    validate_graph_sidecar,
)


def positive_dg0_fixture():
    graph = nx.DiGraph()
    graph.add_node("kegg:C1", group="Compound", name="C1")
    graph.add_node("kegg:C2", group="Compound", name="C2")
    graph.add_node("kegg:R1", group="Reaction", name="R1")
    graph.add_edge("kegg:C2", "kegg:R1")
    graph.add_edge("kegg:R1", "kegg:C1")
    row = {
        "reaction_node": "kegg:R1",
        "mnx_left_compounds": "C1",
        "mnx_right_compounds": "C2",
        "input_compounds": "C2",
        "output_compounds": "C1",
        "input_smiles": "O",
        "output_smiles": "C",
        "I_to_O_string": "O>>C",
        "O_to_I_string": "C>>O",
        "dg0": "1.5",
        "orientation": "reverse",
        "orientation_source": "dg0",
        "thermo_agrees_with_mnx": False,
        "orientation_flip_from_mnx": True,
    }
    return graph, pd.DataFrame([row])


class ValidationTests(unittest.TestCase):
    def test_valid_positive_dg0_orientation_passes(self):
        graph, sidecar = positive_dg0_fixture()
        self.assertEqual(validate_graph_sidecar(graph, sidecar), [])

    def test_double_flip_disagrees_with_full_graph(self):
        full_graph, _ = positive_dg0_fixture()
        final_graph = nx.DiGraph()
        final_graph.add_nodes_from(full_graph.nodes(data=True))
        final_graph.add_edge("kegg:C1", "kegg:R1")
        final_graph.add_edge("kegg:R1", "kegg:C2")
        double_flipped = pd.DataFrame(
            [
                {
                    "reaction_node": "kegg:R1",
                    "mnx_left_compounds": "C2",
                    "mnx_right_compounds": "C1",
                    "input_compounds": "C1",
                    "output_compounds": "C2",
                    "input_smiles": "C",
                    "output_smiles": "O",
                    "I_to_O_string": "C>>O",
                    "O_to_I_string": "O>>C",
                    "dg0": "1.5",
                    "orientation": "reverse",
                    "orientation_source": "dg0",
                    "thermo_agrees_with_mnx": False,
                    "orientation_flip_from_mnx": True,
                }
            ]
        )
        self.assertEqual(validate_graph_sidecar(final_graph, double_flipped), [])
        self.assertTrue(validate_full_graph_sidecar(full_graph, double_flipped))

    def test_embedding_validation_rejects_nonfinite_vector(self):
        _, sidecar = positive_dg0_fixture()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "reaction_embeddings_drfp.parquet"
            pd.DataFrame(
                {
                    "reaction_node": ["kegg:R1"],
                    "embedding_status": ["ok"],
                    "error_message": [""],
                    "embedding": [np.full(2048, np.nan)],
                }
            ).to_parquet(path, index=False)
            result = validate_embedding_table(sidecar, path, "drfp")
        self.assertFalse(result["passes"])
        self.assertEqual(result["finite_rows"], 0)


if __name__ == "__main__":
    unittest.main()
