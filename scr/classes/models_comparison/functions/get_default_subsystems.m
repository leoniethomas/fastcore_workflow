function [idx_pathways,names_pathways] = get_default_subsystems(project, reference_model)
                gly = find(matches(string(project.models.(reference_model).model.subSystems),"Glycolysis/gluconeogenesis"));
                tca = find(matches(string(project.models.(reference_model).model.subSystems),"Citric acid cycle"));
                PPP = find(matches(string(project.models.(reference_model).model.subSystems),"Pentose phosphate pathway"));
                ex = find(matches(string(project.models.(reference_model).model.subSystems),"Exchange/demand reaction"));
                pyr = find(matches(string(project.models.(reference_model).model.subSystems),"Pyruvate metabolism"));
                purine = find(contains(string(project.models.(reference_model).model.subSystems),"Purine "));
            
                pyrimidine =find(contains(string(project.models.(reference_model).model.subSystems),"Pyrimidine "));
            
                nuc = find(matches(string(project.models.(reference_model).model.subSystems),"Nucleotide interconversion"));
                glut = find(matches(string(project.models.(reference_model).model.subSystems),"Glutamate metabolism"));
                Urea_cycle = find(matches(string(project.models.(reference_model).model.subSystems),"Urea cycle"));
                
                proline = find(matches(string(project.models.(reference_model).model.subSystems),"Arginine and proline metabolism"));
                
                % pick amino acids and lipids as one system
                subs = string(project.models.(reference_model).model.subSystems);
                mask = contains(lower(subs), ...
                    ["alanine","glycine","valine","leucine","isoleucine","serine","threonine","cysteine","methionine","aspartate","asparagine","glutamate","glutamine","arginine","proline","histidine","phenylalanine","tyrosine","tryptophan"]);
                AA = find(mask);
                
                lipid_subsystems = [
                    "Fatty acid oxidation"
                    "Fatty acid synthesis"
                    "Glycerophospholipid metabolism"
                    "Sphingolipid metabolism"
                    "Cholesterol metabolism"
                ];
                
                mask = ismember(subs, lipid_subsystems);
                Lipids = find(mask);
            
                idx_pathways = {gly,tca,PPP,ex, pyr,purine, pyrimidine,nuc,glut,Urea_cycle,proline,Lipids};
                names_pathways = ["Glycolysis/gluconeogenesis","Citric acid cycle","Pentose phosphate pathway","Exchange/demand reaction","Pyruvate metabolism",...
                                  "Purine metabolism", "Pyrimidine metabolism", "Nucleotide interconversion", "Glutamate metabolism","Urea cycle",...
                                  "Arginine and proline metabolism","Lipid metabolism"];
end
