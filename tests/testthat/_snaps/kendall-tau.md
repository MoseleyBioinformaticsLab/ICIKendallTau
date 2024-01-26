# large kendall returns correct

    Code
      ici_val
    Output
              tau      pvalue     tau_max 
      -0.00123518  0.67867094  1.00000000 

# completeness works correctly

    Code
      x_comp[4:6, ]
    Output
        s1 s2 core missingness completeness
      4 s1 s5    1           0         1.00
      5 s1 s6    1           2         0.96
      6 s1 s7    1           1         0.98

# kt_fast works properly

    Code
      na_pairs_everything[c("tau", "pvalue")]
    Output
         tau pvalue 
          NA     NA 

---

    Code
      na_pairs_complete[c("tau", "pvalue")]
    Output
              tau      pvalue 
      0.003092146 0.963830701 

---

    Code
      na_matrix_complete[c("tau", "pvalue")]
    Output
      $tau
                  s1           s2           s3         s4
      s1 1.000000000 0.0030921459 0.0072150072 0.10904968
      s2 0.003092146 1.0000000000 0.0006184292 0.04555762
      s3 0.007215007 0.0006184292 1.0000000000 0.01669759
      s4 0.109049680 0.0455576170 0.0166975881 1.00000000
      
      $pvalue
                   s1           s2           s3           s4
      s1 1.076521e-48 9.638307e-01 9.157333e-01 1.097676e-01
      s2 9.638307e-01 1.076521e-48 9.927638e-01 5.040615e-01
      s3 9.157333e-01 9.927638e-01 1.076521e-48 8.065540e-01
      s4 1.097676e-01 5.040615e-01 8.065540e-01 1.076521e-48
      

---

    Code
      na_matrix_pairwise[c("tau", "pvalue")]
    Output
      $tau
                  s1          s2          s3         s4
      s1 1.000000000 0.003092146 0.007215007 0.10904968
      s2 0.003092146 1.000000000 0.002424242 0.04444444
      s3 0.007215007 0.002424242 1.000000000 0.01010101
      s4 0.109049680 0.044444444 0.010101010 1.00000000
      
      $pvalue
                   s1           s2           s3           s4
      s1 1.076521e-48 9.638307e-01 9.157333e-01 1.097676e-01
      s2 9.638307e-01 3.480281e-49 9.714917e-01 5.123482e-01
      s3 9.157333e-01 9.714917e-01 3.480281e-49 8.816279e-01
      s4 1.097676e-01 5.123482e-01 8.816279e-01 3.480281e-49
      

