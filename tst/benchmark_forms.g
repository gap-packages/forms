GOMinus := function(nn, qq)
    return GeneralOrthogonalGroup(-1, nn, qq);
end;

RandomMatrixInvertible := function(n, q)
    local F, M, i, j;
    F := GF(q);
    while true do
        M := NullMat(n, n, F);
        for i in [1..n] do
            for j in [1..n] do
                M[i][j] := PseudoRandom(F);
            od;
        od;
        ConvertToMatrixRep(M, F);
        if RankMat(M) = n then
            return M;
        fi;
    od;
end;

FileBenchmark := function()
    local bench_time, number_to_average, Ns, Qs, Gs, Gs_names, n, q, k, g, G, mine_total, theirs_total, i, H_conj, H, start_mine, time_after, time_start, time_after_after, average_mine, average_theirs, waste1, waste2, file_name;
    bench_time := Runtime();
    number_to_average := 10;
    Ns := [100, 300, 500];
    Qs := [3, 5, 3^2];
    file_name := "benchmark.csv";
    Gs := [Sp, SU, GOMinus];
    Gs_names := ["Sp", "SU", "GOMinus"];
    PrintTo(file_name, "GroupName;n;q;average_time_ours;average_time_theirs");
    AppendTo(file_name, "\n");
    for n in Ns do
        for q in Qs do
            for k in [1..Size(Gs)] do
                g := Gs[k];
                G := g(n, q);
                mine_total := 0;
                theirs_total := 0;
                for i in [1..number_to_average] do
                    H_conj := RandomMatrixInvertible(n, q);
                    H := G^H_conj;
                    start_mine := Runtime();
                    waste1 := PreservedFormsWithScalars(H);
                    time_after := Runtime() - start_mine;
                    time_start := Runtime();
                    waste2 := PreservedForms(H);
                    time_after_after := Runtime() - time_start;
                    mine_total := mine_total + time_after;
                    theirs_total := theirs_total + time_after_after;
                od;
                # '*1000' to convert to seconds
                average_mine := Float(mine_total/(number_to_average*1000));
                average_theirs := Float(theirs_total/(number_to_average*1000));
                AppendTo(file_name, Gs_names[k],";",n,";",q,";",average_mine,";",average_theirs,"\n");
            od;
        od;
    od;
    Print("Benchmark took ", Float((Runtime() - bench_time)/1000), " seconds.");
end;
