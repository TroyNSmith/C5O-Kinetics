-- Reference Tables
CREATE TABLE SMILES (
    smiles_id INTEGER PRIMARY KEY,
    smiles_text TEXT NOT NULL UNIQUE,
    multiplicity INTEGER NOT NULL,
    initial TEXT NOT NULL
);

CREATE TABLE METHODS (
    method_id INTEGER PRIMARY KEY,
    functional TEXT NOT NULL,
    basis TEXT NOT NULL,
    method TEXT NOT NULL,
    inp_template TEXT NOT NULL,
    submit_template TEXT NOT NULL,
    UNIQUE (functional, basis, method)
);

-- Dependent Tables
CREATE TABLE Energies (
    energy_id INTEGER PRIMARY KEY,
    smiles_id INTEGER REFERENCES SMILES(smiles_id),
    method_id INTEGER REFERENCES Methods(method_id),
    energy_value DOUBLE PRECISION NOT NULL,
    units TEXT DEFAULT 'kcal/mol'
);

CREATE TABLE XYZs (
    xyz_id INTEGER PRIMARY KEY,
    smiles_id INTEGER REFERENCES SMILES(smiles_id),
    method_id INTEGER REFERENCES Methods(method_id),
    xyz_text TEXT NOT NULL
);

CREATE TABLE Hessians (
    hessian_id INTEGER PRIMARY KEY,
    smiles_id INTEGER REFERENCES SMILES(smiles_id),
    method_id INTEGER REFERENCES Methods(method_id),
    hessian_data TEXT
);

CREATE TABLE ImaginaryFrequencies (
    freq_id INTEGER PRIMARY KEY,
    smiles_id INTEGER REFERENCES SMILES(smiles_id),
    method_id INTEGER REFERENCES Methods(method_id),
    frequency DOUBLE PRECISION NOT NULL,
    xyz TEXT NOT NULL
);

CREATE TABLE Scans (
    scan_id INTEGER PRIMARY KEY,
    smiles_id INTEGER REFERENCES SMILES(smiles_id),
    method_id INTEGER REFERENCES Methods(method_id),
    scan_data TEXT NOT NULL
);