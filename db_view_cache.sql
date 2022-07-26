DROP PROCEDURE IF EXISTS cache_builder;
DELIMITER ;;

CREATE PROCEDURE cache_builder()
BEGIN
    DROP TABLE IF EXISTS pfamA_pfamseq;
    CREATE VIEW pfamA_pfamseq AS
    SELECT CONCAT(pfamseq.pfamseq_id, '/', CAST(pfamA_reg_full_significant.seq_start AS CHAR), '-', CAST(pfamA_reg_full_significant.seq_end AS CHAR)) AS 'pfamseq_id',
            pfamA_reg_full_significant.pfamA_acc AS 'pfamA_Acc',
            pfamA_reg_full_significant.pfamseq_acc AS 'pfamseq_acc',
            pfamA_reg_full_significant.pfamseq_acc IN (
                SELECT DISTINCT pdb_pfamA_reg.pfamseq_acc
                FROM pdb_pfamA_reg
                WHERE pdb_pfamA_reg.pfamA_acc = pfamA_acc) AS 'has_pdb'
    FROM pfamseq
    LEFT JOIN pfamA_reg_full_significant
        ON pfamA_reg_full_significant.in_full = 1
        AND pfamA_reg_full_significant.pfamseq_acc = pfamseq.pfamseq_acc
    WHERE pfamA_reg_full_significant.pfamA_acc = pfamA_acc;
END;
;;

CALL cache_builder();
