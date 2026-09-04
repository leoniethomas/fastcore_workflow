function pathways = getEssentialPathwayMetabolites(project,reference_model)
    % source for pathway metabolite composition: https://github.com/sysbiolux/MidbrainOrganoid_Miro1_scMetabMod/blob/main/ModelAnalysis/compareFBA_v3_PAPER.m
    
    arguments
        project
        reference_model
    end
        
        pathways.Glycolysis = ["glc_D[c]","g6p[c]","f6p[c]","fdp[c]","dhap[c]","g3p[c]","13dpg[c]","3pg[c]","2pg[c]","pep[c]","pyr[c]"];
        pathways.PPP=["ru5p_D[c]","s7p[c]","e4p[c]"];
        pathways.TCA_cytoplasm_only= ["cit[c]","icit[c]","akg[c]","succ[c]","fum[c]","oaa[c]"];
        pathways.Tricarboxylic_acid_cycle=["cit[m]","icit[m]","akg[m]","succoa[m]","succ[m]","fum[m]","mal_L[m]","oaa[m]"];
        pathways.Oxidative_phosphorylation=["nadh[m]","fadh2[m]","focytC[m]","q10h2[m]","atp[m]"];
        pathways.Fatty_acid_oxidation=["accoa[c]","accoa[m]","coa[c]","coa[m]"];
        %%
        temp=find(ismember(string(project.models.(reference_model).model.subSystems),'Fatty acid oxidation'));
        metList=findMetsFromRxns(project.models.(reference_model).model, ...
                                 project.models.(reference_model).model.rxns(temp))';
        temp=find(contains(metList,'coa'));
        pathways.Fatty_acid_oxidation__all__coa=string(metList(temp));
        %%
        temp=find(ismember(string(project.models.(reference_model).model.subSystems),'Fatty acid oxidation'));
        metList=findMetsFromRxns(project.models.(reference_model).model, ...
                                 project.models.(reference_model).model.rxns(temp))';
        temp=find(contains(metList,'crn'));
        pathways.Fatty_acid_oxidation__all__crn=string(metList(temp));
        %%
        pathways.Cholesterol=["sql[r]","zymst[r]","chsterol[r]","chsterol[c]"];
        pathways.Dopamine_Tyrosine_metabolism=["phe_L[c]","tyr_L[c]","34dhphe[c]","dopa[c]","34dhpac[c]","34dhpha[c]","34dhpe[c]","tym[c]"];

        essential_pathway_metabolites = pathways;
end