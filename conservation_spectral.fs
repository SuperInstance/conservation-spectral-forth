\ ============================================================
\ conservation_spectral.fs — Conservation Spectral SDK in Forth
\ Stack-Based = Sheaf-Theoretic Implementation
\
\ The deep insight: Forth's stack is the STALK. Each word is a
\ LOCAL SECTION. Composition via concatenation is RESTRICTION.
\ FULL-ANALYSIS is the GLOBAL SECTION gluing all local sections.
\
\ Compatible with Gforth. Version 0.1.0
\ ============================================================

10 CONSTANT MAX-N
MAX-N DUP * CONSTANT MAX-MATRIX

CREATE TRANSITION MAX-MATRIX FLOATS ALLOT
CREATE LAPLACIAN  MAX-MATRIX FLOATS ALLOT
CREATE EIGENVALS  MAX-N FLOATS ALLOT
CREATE EIGENVECS  MAX-MATRIX FLOATS ALLOT
CREATE WORK-VEC   MAX-N FLOATS ALLOT
CREATE WORK-VEC2  MAX-N FLOATS ALLOT
CREATE ATTRIBUTES MAX-N FLOATS ALLOT
CREATE DEFLATED   MAX-MATRIX FLOATS ALLOT

VARIABLE GRAPH-SIZE

\ Temp float storage (4 slots)
CREATE F-TEMP 4 FLOATS ALLOT
: F-TEMP@ FLOATS F-TEMP F@ ;
: F-TEMP! FLOATS F-TEMP F! ;

\ === Matrix Operations ===
\ Row-major: index = row * n + col
\ Stack: ( row col base -- addr )
: MATRIX-ADDR ( row col base -- addr )
    ROT GRAPH-SIZE @ * ROT + FLOATS + ;

: MATRIX-GET ( row col base -- f:value )
    MATRIX-ADDR F@ ;

: MATRIX-SET ( f:value row col base -- )
    MATRIX-ADDR F! ;

\ Clear n×n matrix
: CLEAR-MATRIX ( base -- )
    GRAPH-SIZE @ DUP * 0 DO
        0E0 OVER I FLOATS + F!
    LOOP
    DROP ;

\ === Vector operations on WORK-VEC ===
: VEC-NORM ( -- f:norm )
    0E0
    GRAPH-SIZE @ 0 DO
        WORK-VEC I FLOATS + F@ FDUP F* F+
    LOOP
    FSQRT ;

: VEC-NORMALIZE ( -- )
    VEC-NORM
    FDUP 0E0 F> IF
        GRAPH-SIZE @ 0 DO
            WORK-VEC I FLOATS + F@ FOVER F/
            WORK-VEC I FLOATS + F!
        LOOP
    THEN
    FDROP ;

: VEC-COPY2TO1 ( -- )
    GRAPH-SIZE @ 0 DO
        WORK-VEC2 I FLOATS + F@
        WORK-VEC I FLOATS + F!
    LOOP ;

\ Matrix-vector multiply: WORK-VEC2 = Matrix × WORK-VEC
\ base is the matrix address
VARIABLE MATVEC-ROW
: MATVEC ( base -- )
    GRAPH-SIZE @ 0 DO
        I MATVEC-ROW !
        0E0
        GRAPH-SIZE @ 0 DO
            WORK-VEC I FLOATS + F@
            MATVEC-ROW @ I OVER MATRIX-GET
            F* F+
        LOOP
        WORK-VEC2 MATVEC-ROW @ FLOATS + F!
    LOOP
    DROP ;

VARIABLE LAP-ROW

\ ============================================================
\ 1. BUILD-LAPLACIAN
\    L = D - W from TRANSITION adjacency matrix
\ ============================================================
: BUILD-LAPLACIAN ( -- )
    LAPLACIAN CLEAR-MATRIX
    GRAPH-SIZE @ 0 DO
        I LAP-ROW !
        \ Compute degree of vertex I (row sum)
        0E0
        GRAPH-SIZE @ 0 DO
            LAP-ROW @ I TRANSITION MATRIX-GET F+
        LOOP
        \ FP stack: degree
        GRAPH-SIZE @ 0 DO
            LAP-ROW @ I = IF
                FDUP
                LAP-ROW @ I TRANSITION MATRIX-GET F-
                LAP-ROW @ I LAPLACIAN MATRIX-SET
            ELSE
                LAP-ROW @ I TRANSITION MATRIX-GET FNEGATE
                LAP-ROW @ I LAPLACIAN MATRIX-SET
            THEN
        LOOP
        FDROP
    LOOP ;

\ ============================================================
\ 2. POWER-ITERATION
\    Finds largest eigenvalue of matrix at base
\    Output: eigenvalue on FP stack, eigenvector in WORK-VEC
\ ============================================================
: POWER-ITERATION ( base -- f:eigenvalue )
    { base }
    GRAPH-SIZE @ 0 DO
        1E0 WORK-VEC I FLOATS + F!
    LOOP
    VEC-NORMALIZE

    0E0 \ previous eigenvalue
    500 0 DO
        base MATVEC
        VEC-COPY2TO1
        VEC-NORMALIZE
        base MATVEC
        0E0
        GRAPH-SIZE @ 0 DO
            WORK-VEC I FLOATS + F@
            WORK-VEC2 I FLOATS + F@
            F* F+
        LOOP
        FOVER F- FABS 1E-9 F< IF
            NIP UNLOOP EXIT
        THEN
    LOOP
    NIP ;

\ ============================================================
\ 3. INVERSE-POWER — finds smallest eigenvalue of LAPLACIAN
\ ============================================================
: INVERSE-POWER ( -- f:smallest_eigenvalue )
    0E0
    GRAPH-SIZE @ 0 DO
        LAPLACIAN I I MATRIX-GET FMAX
    LOOP
    0 F-TEMP!  \ save shift

    GRAPH-SIZE @ 0 DO
        GRAPH-SIZE @ 0 DO
            0 F-TEMP@ LAPLACIAN I J MATRIX-GET F-
            I J = IF 0 F-TEMP@ F+ THEN
            DEFLATED I J MATRIX-SET
        LOOP
    LOOP

    DEFLATED POWER-ITERATION
    0 F-TEMP@ FSWAP F- ;

\ ============================================================
\ COMPUTE-EIGENVALUES (all eigenvalues via deflation)
\ ============================================================
: COMPUTE-EIGENVALUES ( -- )
    0E0
    GRAPH-SIZE @ 0 DO
        LAPLACIAN I I MATRIX-GET FMAX
    LOOP
    0 F-TEMP!  \ shift

    \ Build M = shift*I - L in DEFLATED
    GRAPH-SIZE @ 0 DO
        GRAPH-SIZE @ 0 DO
            0 F-TEMP@ LAPLACIAN I J MATRIX-GET F-
            I J = IF 0 F-TEMP@ F+ THEN
            DEFLATED I J MATRIX-SET
        LOOP
    LOOP

    GRAPH-SIZE @ 0 DO
        \ Seed vector
        GRAPH-SIZE @ 0 DO
            I 1+ S>F WORK-VEC I FLOATS + F!
        LOOP
        VEC-NORMALIZE

        0E0
        300 0 DO
            DEFLATED MATVEC
            VEC-COPY2TO1
            VEC-NORMALIZE
            DEFLATED MATVEC
            0E0
            GRAPH-SIZE @ 0 DO
                WORK-VEC I FLOATS + F@
                WORK-VEC2 I FLOATS + F@
                F* F+
            LOOP
            FOVER F- FABS 1E-8 F< IF
                NIP UNLOOP
            THEN
        LOOP
        NIP

        \ eigenvalue of L = shift - eigenvalue_M
        0 F-TEMP@ FSWAP F-
        EIGENVALS I FLOATS + F!

        \ Deflate
        GRAPH-SIZE @ 0 DO
            GRAPH-SIZE @ 0 DO
                WORK-VEC I FLOATS + F@
                WORK-VEC J FLOATS + F@ F*
                DEFLATED I J DUP MATRIX-GET F-
                DEFLATED I J MATRIX-SET
            LOOP
        LOOP
    LOOP ;

\ Sort eigenvalues ascending
: SORT-EIGENVALUES ( -- )
    GRAPH-SIZE @ 1- 0 DO
        GRAPH-SIZE @ 1- I ?DO
            EIGENVALS J FLOATS + F@
            EIGENVALS J 1+ FLOATS + F@
            F> IF
                EIGENVALS J FLOATS + F@
                EIGENVALS J 1+ FLOATS + F@
                EIGENVALS J FLOATS + F!
                EIGENVALS J 1+ FLOATS + F!
            THEN
        LOOP
    LOOP ;

\ ============================================================
\ 4. CONSERVATION-RATIO
\    Variance of gradient of projected attributes
\ ============================================================
: CONSERVATION-RATIO ( k -- f:ratio )
    DROP
    GRAPH-SIZE @ 1- 0 DO
        ATTRIBUTES I FLOATS + F@
        ATTRIBUTES I 1+ FLOATS + F@ F-
        WORK-VEC I FLOATS + F!
    LOOP

    0E0
    GRAPH-SIZE @ 1- 0 DO
        WORK-VEC I FLOATS + F@ F+
    LOOP
    GRAPH-SIZE @ 1- S>F F/
    1 F-TEMP!

    0E0
    GRAPH-SIZE @ 1- 0 DO
        WORK-VEC I FLOATS + F@ 1 F-TEMP@ F- FDUP F* F+
    LOOP
    GRAPH-SIZE @ 1- S>F F/ ;

\ ============================================================
\ 5. SPECTRAL-GAP
\    Largest gap between consecutive eigenvalues
\ ============================================================
: SPECTRAL-GAP ( -- f:gap )
    0E0
    GRAPH-SIZE @ 1- 0 DO
        EIGENVALS I 1+ FLOATS + F@
        EIGENVALS I FLOATS + F@
        F- FMAX
    LOOP ;

\ ============================================================
\ 6. CHEEGER-APPROX
\    Approximate Cheeger constant from Fiedler vector sign cut
\ ============================================================
: CHEEGER-APPROX ( -- f:cheeger )
    0E0 2 F-TEMP!  \ cut
    0E0 3 F-TEMP!  \ vol_s

    GRAPH-SIZE @ 0 DO
        WORK-VEC I FLOATS + F@ F0< IF
            LAPLACIAN I I MATRIX-GET 3 F-TEMP@ F+ 3 F-TEMP!
            GRAPH-SIZE @ 0 DO
                I J <> IF
                    WORK-VEC J FLOATS + F@ F0>= IF
                        I J TRANSITION MATRIX-GET 2 F-TEMP@ F+ 2 F-TEMP!
                    THEN
                THEN
            LOOP
        THEN
    LOOP

    0E0
    GRAPH-SIZE @ 0 DO
        LAPLACIAN I I MATRIX-GET F+
    LOOP
    1 F-TEMP!  \ total_vol

    1 F-TEMP@ 3 F-TEMP@ F-  \ vol_comp
    3 F-TEMP@ FMIN
    FDUP 1E-10 F< IF
        FDROP 0E0
    ELSE
        2 F-TEMP@ FSWAP F/
    THEN ;

\ ============================================================
\ 7. ANOMALY-DETECT
\    0=nominal, 1=warning, 2=critical
\ ============================================================
: ANOMALY-DETECT ( f:obs f:mean f:std -- status )
    1 F-TEMP!  \ std
    FDUP 1E-15 F< IF
        2DROP DROP 0 EXIT
    THEN
    F-  \ obs - mean
    FABS
    1 F-TEMP@ F/
    FDUP 3E0 F> IF FDROP 2 EXIT THEN
    2E0 F> IF DROP 1 EXIT THEN
    DROP 0 ;

\ ============================================================
\ Sliding window tracker
\ ============================================================
10 CONSTANT TRACKER-WINDOW
CREATE TRACKER-HISTORY TRACKER-WINDOW FLOATS ALLOT
VARIABLE TRACKER-COUNT
VARIABLE TRACKER-BASELINE-SET
CREATE TRACKER-MEAN 1 FLOATS ALLOT
CREATE TRACKER-STD  1 FLOATS ALLOT

: TRACKER-RESET ( -- )
    0 TRACKER-COUNT !
    0 TRACKER-BASELINE-SET !
    TRACKER-WINDOW 0 DO
        0E0 TRACKER-HISTORY I FLOATS + F!
    LOOP ;

: TRACKER-FEED ( f:observation -- status )
    TRACKER-COUNT @ TRACKER-WINDOW < IF
        TRACKER-HISTORY TRACKER-COUNT @ FLOATS + F!
        1 TRACKER-COUNT +!
        TRACKER-COUNT @ TRACKER-WINDOW = IF
            0E0
            TRACKER-WINDOW 0 DO
                TRACKER-HISTORY I FLOATS + F@ F+
            LOOP
            TRACKER-WINDOW S>F F/
            TRACKER-MEAN F!
            0E0
            TRACKER-WINDOW 0 DO
                TRACKER-HISTORY I FLOATS + F@
                TRACKER-MEAN F@ F- FDUP F* F+
            LOOP
            TRACKER-WINDOW S>F F/
            FSQRT TRACKER-STD F!
            1 TRACKER-BASELINE-SET !
            0
        ELSE
            0
        THEN
    ELSE
        TRACKER-WINDOW 1- 0 DO
            TRACKER-HISTORY I 1+ FLOATS + F@
            TRACKER-HISTORY I FLOATS + F!
        LOOP
        TRACKER-HISTORY TRACKER-WINDOW 1- FLOATS + F!
        TRACKER-HISTORY TRACKER-WINDOW 1- FLOATS + F@
        TRACKER-MEAN F@
        TRACKER-STD F@
        ANOMALY-DETECT
    THEN ;

\ ============================================================
\ 8. FULL-ANALYSIS — The GLOBAL SECTION
\ ============================================================
: FULL-ANALYSIS ( -- )
    CR ." === Conservation Spectral Analysis ===" CR
    ." Graph size: " GRAPH-SIZE @ . CR

    BUILD-LAPLACIAN
    ." Laplacian built." CR

    COMPUTE-EIGENVALUES
    SORT-EIGENVALUES
    ." Eigenvalues:" CR
    GRAPH-SIZE @ 0 DO
        ."   λ[" I . ." ] = "
        EIGENVALS I FLOATS + F@ F. CR
    LOOP

    ." Spectral gap: " SPECTRAL-GAP F. CR

    GRAPH-SIZE @ 1 > IF
        ." Fiedler value: " EIGENVALS 1 FLOATS + F@ F. CR
        ." Conservation ratio: " 1 CONSERVATION-RATIO F. CR
    THEN

    ." === Analysis Complete ===" CR ;

\ ============================================================
\ Demo graphs
\ ============================================================
: DEMO-GRAPH ( -- )
    4 GRAPH-SIZE !
    TRANSITION CLEAR-MATRIX
    1E0 0 1 TRANSITION MATRIX-SET  1E0 1 0 TRANSITION MATRIX-SET
    1E0 1 2 TRANSITION MATRIX-SET  1E0 2 1 TRANSITION MATRIX-SET
    1E0 2 3 TRANSITION MATRIX-SET  1E0 3 2 TRANSITION MATRIX-SET
    1E0 ATTRIBUTES 0 FLOATS + F!
    2E0 ATTRIBUTES 1 FLOATS + F!
    3E0 ATTRIBUTES 2 FLOATS + F!
    4E0 ATTRIBUTES 3 FLOATS + F!
    ." 4-vertex path graph loaded." CR ;

: DEMO-CYCLE ( -- )
    5 GRAPH-SIZE !
    TRANSITION CLEAR-MATRIX
    1E0 0 1 TRANSITION MATRIX-SET  1E0 1 0 TRANSITION MATRIX-SET
    1E0 1 2 TRANSITION MATRIX-SET  1E0 2 1 TRANSITION MATRIX-SET
    1E0 2 3 TRANSITION MATRIX-SET  1E0 3 2 TRANSITION MATRIX-SET
    1E0 3 4 TRANSITION MATRIX-SET  1E0 4 3 TRANSITION MATRIX-SET
    1E0 4 0 TRANSITION MATRIX-SET  1E0 0 4 TRANSITION MATRIX-SET
    1E0 ATTRIBUTES 0 FLOATS + F!
    2E0 ATTRIBUTES 1 FLOATS + F!
    3E0 ATTRIBUTES 2 FLOATS + F!
    5E0 ATTRIBUTES 3 FLOATS + F!
    8E0 ATTRIBUTES 4 FLOATS + F!
    ." 5-vertex cycle graph loaded." CR ;

\ ============================================================
\ Self-test suite
\ ============================================================
VARIABLE TESTS-PASSED
VARIABLE TESTS-FAILED
: TEST-PASS 1 TESTS-PASSED +! ;
: TEST-FAIL 1 TESTS-FAILED +! ;
: TEST IF TEST-PASS ELSE TEST-FAIL THEN ;

: RUN-TESTS ( -- )
    0 TESTS-PASSED !
    0 TESTS-FAILED !
    CR ." --- Running Self Tests ---" CR

    \ Test 1: Graph size
    4 GRAPH-SIZE !
    GRAPH-SIZE @ 4 = TEST
    ." Test 1 (graph size): " GRAPH-SIZE @ 4 = IF ." PASS" ELSE ." FAIL" THEN CR

    \ Test 2: Matrix set/get
    TRANSITION CLEAR-MATRIX
    3.14E0 1 2 TRANSITION MATRIX-SET
    1 2 TRANSITION MATRIX-GET 3.14E0 0.001E0 F~ TEST
    ." Test 2 (matrix set/get): " 1 2 TRANSITION MATRIX-GET 3.14E0 0.001E0 F~ IF ." PASS" ELSE ." FAIL" THEN CR

    \ Test 3: Laplacian for path graph
    DEMO-GRAPH
    BUILD-LAPLACIAN
    LAPLACIAN 0 0 MATRIX-GET 1E0 0.001E0 F~ TEST
    ." Test 3 (L[0][0]=1): " LAPLACIAN 0 0 MATRIX-GET 1E0 0.001E0 F~ IF ." PASS" ELSE ." FAIL" THEN CR
    LAPLACIAN 1 1 MATRIX-GET 2E0 0.001E0 F~ TEST
    ." Test 3 (L[1][1]=2): " LAPLACIAN 1 1 MATRIX-GET 2E0 0.001E0 F~ IF ." PASS" ELSE ." FAIL" THEN CR
    LAPLACIAN 0 1 MATRIX-GET -1E0 0.001E0 F~ TEST
    ." Test 3 (L[0][1]=-1): " LAPLACIAN 0 1 MATRIX-GET -1E0 0.001E0 F~ IF ." PASS" ELSE ." FAIL" THEN CR
    LAPLACIAN 0 2 MATRIX-GET FABS 0.001E0 F< TEST
    ." Test 3 (L[0][2]=0): " LAPLACIAN 0 2 MATRIX-GET FABS 0.001E0 F< IF ." PASS" ELSE ." FAIL" THEN CR

    \ Test 4: Eigenvalues
    COMPUTE-EIGENVALUES
    SORT-EIGENVALUES
    ." Eigenvalues: " CR
    GRAPH-SIZE @ 0 DO
        ."   λ[" I . ." ] = " EIGENVALS I FLOATS + F@ F. CR
    LOOP
    EIGENVALS 0 FLOATS + F@ FABS 0.1E0 F< TEST
    ." Test 4 (λ[0]≈0): " EIGENVALS 0 FLOATS + F@ FABS 0.1E0 F< IF ." PASS" ELSE ." FAIL" THEN CR
    EIGENVALS 1 FLOATS + F@ FABS 0.3E0 F> TEST
    ." Test 4 (λ[1]>0.3): " EIGENVALS 1 FLOATS + F@ FABS 0.3E0 F> IF ." PASS" ELSE ." FAIL" THEN CR

    \ Test 5: Spectral gap > 0
    SPECTRAL-GAP 0E0 F> TEST
    ." Test 5 (gap>0): " SPECTRAL-GAP 0E0 F> IF ." PASS" ELSE ." FAIL" THEN CR

    \ Test 6: Anomaly detection
    1.5E0 1.0E0 0.5E0 ANOMALY-DETECT 0 = TEST
    ." Test 6 (nominal z=1): " 1.5E0 1.0E0 0.5E0 ANOMALY-DETECT 0 = IF ." PASS" ELSE ." FAIL" THEN CR
    2.5E0 1.0E0 0.5E0 ANOMALY-DETECT 1 = TEST
    ." Test 6 (warning z=3): " 2.5E0 1.0E0 0.5E0 ANOMALY-DETECT 1 = IF ." PASS" ELSE ." FAIL" THEN CR
    5.0E0 1.0E0 0.5E0 ANOMALY-DETECT 2 = TEST
    ." Test 6 (critical z=8): " 5.0E0 1.0E0 0.5E0 ANOMALY-DETECT 2 = IF ." PASS" ELSE ." FAIL" THEN CR

    \ Test 7: Tracker
    TRACKER-RESET
    1.0E0 TRACKER-FEED DROP
    1.1E0 TRACKER-FEED DROP
    0.9E0 TRACKER-FEED DROP
    1.0E0 TRACKER-FEED DROP
    1.05E0 TRACKER-FEED DROP
    0.95E0 TRACKER-FEED DROP
    1.0E0 TRACKER-FEED DROP
    1.1E0 TRACKER-FEED DROP
    0.9E0 TRACKER-FEED DROP
    1.0E0 TRACKER-FEED DROP
    TRACKER-BASELINE-SET @ 1 = TEST
    ." Test 7 (baseline set): " TRACKER-BASELINE-SET @ 1 = IF ." PASS" ELSE ." FAIL" THEN CR
    1.0E0 TRACKER-FEED 0 = TEST
    ." Test 7 (nominal): " 1.0E0 TRACKER-FEED 0 = IF ." PASS" ELSE ." FAIL" THEN CR

    CR ." Passed: " TESTS-PASSED @ . ."  Failed: " TESTS-FAILED @ . CR
    ." --- Tests Complete ---" CR ;

\ ============================================================
\ Main
\ ============================================================
: MAIN ( -- )
    CR ." ╔═══════════════════════════════════════════════════════╗" CR
    ." ║ Conservation Spectral SDK — Forth (Sheaf Edition)    ║" CR
    ." ║ Stack = Stalk | Words = Sections | Glue = Compose    ║" CR
    ." ╚═══════════════════════════════════════════════════════╝" CR

    RUN-TESTS

    CR ." --- Demo: Path Graph ---" CR
    DEMO-GRAPH
    FULL-ANALYSIS

    CR ." --- Demo: Cycle Graph ---" CR
    DEMO-CYCLE
    FULL-ANALYSIS ;

MAIN
BYE
