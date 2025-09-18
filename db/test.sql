CREATE TABLE IF NOT EXISTS `PR` (
	`pr_id` integer primary key NOT NULL UNIQUE,
	`smiles` TEXT NOT NULL UNIQUE,
	`charge` INTEGER NOT NULL,
	`multiplicity` INTEGER NOT NULL,
	`nodes` TEXT NOT NULL,
	`adjacency` TEXT NOT NULL,
	`goat_method` INTEGER NOT NULL,
	`opt_numfreq_method` INTEGER NOT NULL,
	`calc_method` INTEGER NOT NULL,
FOREIGN KEY(`goat_method`) REFERENCES `METHODS`(`method_id`),
FOREIGN KEY(`opt_numfreq_method`) REFERENCES `METHODS`(`method_id`),
FOREIGN KEY(`calc_method`) REFERENCES `METHODS`(`method_id`)
);
CREATE TABLE IF NOT EXISTS `TS` (
	`ts_id` integer primary key NOT NULL UNIQUE,
	`pr_1` INTEGER NOT NULL,
	`pr_2` INTEGER NOT NULL,
	`matrix` TEXT NOT NULL,
	`mechanism` TEXT NOT NULL,
FOREIGN KEY(`pr_1`) REFERENCES `PR`(`pr_id`),
FOREIGN KEY(`pr_2`) REFERENCES `PR`(`pr_id`)
);
CREATE TABLE IF NOT EXISTS `GOAT` (
	`pr_id` INTEGER NOT NULL UNIQUE,
	`method_id` INTEGER NOT NULL,
	`global_minimum` TEXT NOT NULL,
FOREIGN KEY(`pr_id`) REFERENCES `PR`(`pr_id`),
FOREIGN KEY(`method_id`) REFERENCES `PR`(`goat_method`)
);
CREATE TABLE IF NOT EXISTS `METHODS` (
	`method_id` integer primary key NOT NULL UNIQUE,
	`functional` TEXT NOT NULL,
	`basis` TEXT
);
CREATE TABLE IF NOT EXISTS `OPT_NUMFREQ` (
	`pr_id` INTEGER NOT NULL UNIQUE,
	`method_id` INTEGER NOT NULL,
	`data` INTEGER NOT NULL,
FOREIGN KEY(`pr_id`) REFERENCES `PR`(`pr_id`),
FOREIGN KEY(`method_id`) REFERENCES `PR`(`opt_numfreq_method`)
);