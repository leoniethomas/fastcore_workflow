function saveStructToFolder(S,outdir)



save_table_with_rownames(S.metadata, fullfile(outdir, 'metadata.csv'));


save_table_with_rownames(S.raw_counts, fullfile(outdir, 'raw_counts.csv'));


save_table_with_rownames(S.norm_counts, fullfile(outdir, 'norm_counts.csv'));


save_table_with_rownames(S.mapping_exp_2_rxns, fullfile(outdir, 'mapping_exp_2_rxns.csv'));


save_table_with_rownames(S.dico, fullfile(outdir, 'dico.csv'));
    
    function save_table_with_rownames(T, fname)
        if ~isempty(T.Properties.RowNames)
            % Add row names as the first column
            T2 = addvars(T, T.Properties.RowNames, 'Before', 1, 'NewVariableNames', 'RowName');
            writetable(T2, fname);
        else
            writetable(T, fname);
        end
    end

end