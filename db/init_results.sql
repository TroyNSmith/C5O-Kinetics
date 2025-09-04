CREATE TABLE calculations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    smiles TEXT,
    method TEXT,
    proced TEXT,
    mult REAL,
    energy REAL,
    xyz TEXT,
    hess TEXT,
    UNIQUE(smiles, method, proced)
);