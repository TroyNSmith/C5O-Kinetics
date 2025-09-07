"""Dictionary of standard bond lengths"""

# Dictionary of A-B single bond lengths
LEN_DCT = {
    ("H", "H"): 0.74,
    ("H", "C"): 1.09,
    ("H", "N"): 1.01,
    ("H", "O"): 0.95,
    ("H", "Cl"): 1.275,
    ("C", "C"): 1.54,
    ("C", "N"): 1.47,
    ("C", "O"): 1.43,
    ("N", "N"): 1.45,
    ("N", "O"): 1.45,
    ("O", "O"): 1.40,
    ("C", "Cl"): 1.74,
    ("Cl", "Cl"): 2.0,
}


def get_ccsdt_parameters(num_heavy_atoms: int) -> dict:
    num_cpus = [0, 4, 4, 12, 12, 16, 20, 20, 20]
    mem_per_cpu = [0, 1000, 1000, 2000, 4000, 4000, 4000, 6000, 12000]
    lscratch_size = [0, 50, 50, 100, 200, 200, 200, 200, 200]
    time_requested = [
        "00:00:00",
        "00:30:00",
        "00:30:00",
        "00:30:00",
        "00:60:00",
        "02:00:00",
        "02:00:00",
        "02:00:00",
        "06:00:00",
    ]
    idx = min(num_heavy_atoms, len(num_cpus) - 1)
    return {
        "[num_cpus]": num_cpus[idx],
        "[mem_per_cpu]": mem_per_cpu[idx],
        "[input_mem_per_cpu]": int(mem_per_cpu[idx] * .75),
        "[lscratch_size]": lscratch_size[idx],
        "[time_requested]": time_requested[idx],
    }
