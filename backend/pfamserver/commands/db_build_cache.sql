DROP PROCEDURE IF EXISTS cache_builder;
DELIMITER ;;

CREATE PROCEDURE cache_builder()
BEGIN
    DECLARE n INT DEFAULT 0;
    DECLARE i INT DEFAULT 0;
    DECLARE pfamA_acc VARCHAR(7) DEFAULT '';

    DROP TABLE IF EXISTS pfamA_pfamseq;
    CREATE TABLE pfamA_pfamseq (
        pfamseq_id VARCHAR(40),
        pfamA_acc VARCHAR(7),
        pfamseq_acc VARCHAR(10),
        has_pdb BOOLEAN NOT NULL DEFAULT 0,
        CONSTRAINT pk_constraint PRIMARY KEY (pfamseq_id, pfamA_acc, pfamseq_acc)
    );

    SELECT COUNT(*) FROM pfamA INTO n;
    SET i=0;
    WHILE i<n DO
        SELECT pfamA.pfamA_acc FROM pfamA ORDER BY pfamA.pfamA_acc LIMIT i,1 INTO pfamA_acc;
        SELECT i, pfamA_acc;
        INSERT INTO pfamA_pfamseq
        SELECT CONCAT(pfamseq.pfamseq_id, '/', CAST(pfamA_reg_full_significant.seq_start AS CHAR), '-', CAST(pfamA_reg_full_significant.seq_end AS CHAR)) AS 'pfamseq_id',
               pfamA_reg_full_significant.pfamA_acc,
               pfamA_reg_full_significant.pfamseq_acc,
               pfamA_reg_full_significant.pfamseq_acc IN (
                    SELECT DISTINCT pdb_pfamA_reg.pfamseq_acc
                    FROM pdb_pfamA_reg
                    WHERE pdb_pfamA_reg.pfamA_acc = pfamA_acc) AS 'has_pdb'
        FROM pfamseq
        LEFT JOIN pfamA_reg_full_significant
            ON pfamA_reg_full_significant.in_full = 1
            AND pfamA_reg_full_significant.pfamseq_acc = pfamseq.pfamseq_acc
        WHERE pfamA_reg_full_significant.pfamA_acc = pfamA_acc;
        SET i = i + 1;
    END WHILE;
    CREATE INDEX pfamA_pfamseq_pfamseq_id_IDX USING HASH ON pfamA_pfamseq (pfamseq_id);
    CREATE INDEX pfamA_pfamseq_pfamA_acc_IDX USING HASH ON pfamA_pfamseq (pfamA_acc);
    CREATE INDEX pfamA_pfamseq_pfamseq_acc_IDX USING HASH ON pfamA_pfamseq (pfamseq_acc);
END;
;;

CALL cache_builder();