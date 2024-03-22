from common import mol_h2o_3, specifications
from utils import generate_component_data

# Generates test data from component
# This does not generate reference data, but the outputs of the individual calculations
# that would go into a QCManybody analysis. This is so you don't need to run full
# computations to test

#################################################
# Single-level, without supersystem
#################################################
generate_component_data(
    mol_h2o_3,
    {1: "e_scf", 2: "e_scf", 3: "e_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_single_energy_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "e_scf", 2: "e_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_single_energy_2",
)

generate_component_data(
    mol_h2o_3,
    {1: "g_scf", 2: "g_scf", 3: "g_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_single_gradient_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "g_scf", 2: "g_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_single_gradient_2",
)

generate_component_data(
    mol_h2o_3,
    {1: "h_scf", 2: "h_scf", 3: "h_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_single_hessian_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "h_scf", 2: "h_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_single_hessian_2",
)

#################################################
# Multi-level, without supersystem
#################################################

generate_component_data(
    mol_h2o_3,
    {1: "e_mp2", 2: "e_b3lyp", 3: "e_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_energy_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "e_mp2", 2: "e_b3lyp"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_energy_2",
)

generate_component_data(
    mol_h2o_3,
    {1: "g_mp2", 2: "g_b3lyp", 3: "g_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_gradient_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "g_mp2", 2: "g_b3lyp"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_gradient_2",
)

generate_component_data(
    mol_h2o_3,
    {1: "h_scf_atz", 2: "h_scf_adz", 3: "h_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_hessian_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "h_scf_atz", 2: "h_scf_adz"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_hessian_2",
)

#################################################
# Multi-level, with supersystem
#################################################
generate_component_data(
    mol_h2o_3,
    {1: "e_mp2", 2: "e_b3lyp", "supersystem": "e_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_ss_energy_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "e_mp2", "supersystem": "e_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_ss_energy_2",
)

generate_component_data(
    mol_h2o_3,
    {1: "g_mp2", 2: "g_b3lyp", "supersystem": "g_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_ss_gradient_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "g_mp2", "supersystem": "g_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_ss_gradient_2",
)

generate_component_data(
    mol_h2o_3,
    {1: "h_scf_atz", 2: "h_scf_adz", "supersystem": "h_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_ss_hessian_1",
)

generate_component_data(
    mol_h2o_3,
    {1: "h_scf_atz", "supersystem": "h_scf"},
    specifications,
    ["cp", "nocp", "vmfc"],
    True,
    "h2o_trimer_multi_ss_hessian_2",
)
