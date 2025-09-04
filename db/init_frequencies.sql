CREATE TABLE calculations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    smiles TEXT,
    method TEXT,
    proced TEXT,
    mult REAL,
    mode TEXT,
    freq TEXT,
    xyz TEXT,
    UNIQUE(smiles, method, proced)
);