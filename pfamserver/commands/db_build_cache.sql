DROP TABLE IF EXISTS pfamA_pfamseq;

CREATE TABLE pfamA_pfamseq ( pfamseq_id VARCHAR(16), pfamA_acc VARCHAR(7), pfamseq_acc VARCHAR(10), has_pdb BOOLEAN NOT NULL DEFAULT 0)
SELECT CONCAT(pfamseq.pfamseq_id, '/', CONVERT(pfamA_reg_full_significant.seq_start, CHAR), '-', CONVERT(pfamA_reg_full_significant.seq_end, CHAR)) AS 'pfamseq_id',
       pfamA_reg_full_significant.pfamA_acc,
       pfamA_reg_full_significant.pfamseq_acc,
       COUNT(pdb_pfamA_reg.pdb_id) AS 'has_pdb'
FROM pfamseq
LEFT JOIN pfamA_reg_full_significant
    ON pfamA_reg_full_significant.in_full = 1
    AND pfamA_reg_full_significant.pfamseq_acc = pfamseq.pfamseq_acc
LEFT JOIN pdb_pfamA_reg
    ON pdb_pfamA_reg.pfamseq_acc = pfamA_reg_full_significant.pfamseq_acc
    AND pdb_pfamA_reg.pfamA_acc = pfamA_reg_full_significant.pfamA_acc
GROUP BY pfamseq.pfamseq_id, pfamA_reg_full_significant.seq_start, pfamA_reg_full_significant.seq_end, pfamA_reg_full_significant.pfamA_acc, pfamA_reg_full_significant.pfamseq_acc;