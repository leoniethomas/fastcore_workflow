function [idxPathways, namesPathways] = getDefaultSubsystems(project, referenceModel)
                gly = find(matches(string(project.models.(referenceModel).model.subSystems), "Glycolysis/gluconeogenesis"));
                tca = find(matches(string(project.models.(referenceModel).model.subSystems), "Citric acid cycle"));
                PPP = find(matches(string(project.models.(referenceModel).model.subSystems), "Pentose phosphate pathway"));
                ex = find(matches(string(project.models.(referenceModel).model.subSystems), "Exchange/demand reaction"));
                pyr = find(matches(string(project.models.(referenceModel).model.subSystems), "Pyruvate metabolism"));
                purine = find(contains(string(project.models.(referenceModel).model.subSystems), "Purine "));
            
                pyrimidine = find(contains(string(project.models.(referenceModel).model.subSystems), "Pyrimidine "));
            
                nuc = find(matches(string(project.models.(referenceModel).model.subSystems), "Nucleotide interconversion"));
                glut = find(matches(string(project.models.(referenceModel).model.subSystems), "Glutamate metabolism"));
                ureaCycle = find(matches(string(project.models.(referenceModel).model.subSystems), "Urea cycle"));
                
                proline = find(matches(string(project.models.(referenceModel).model.subSystems), "Arginine and proline metabolism"));
                
                % pick amino acids and lipids as one system
                subs = string(project.models.(referenceModel).model.subSystems);
                mask = contains(lower(subs), ...
                    ["alanine", "glycine", "valine", "leucine", "isoleucine", "serine", "threonine", "cysteine", "methionine", ...
                    "aspartate", "asparagine", "glutamate", "glutamine", "arginine", "proline", "histidine", "phenylalanine", "tyrosine", "tryptophan"]);
                AA = find(mask); % not used anymore ???
                
                lipidSubsystems = [
                    "Fatty acid oxidation"
                    "Fatty acid synthesis"
                    "Glycerophospholipid metabolism"
                    "Sphingolipid metabolism"
                    "Cholesterol metabolism"
                ];
                
                mask = ismember(subs, lipidSubsystems);
                lipids = find(mask);
            
                idxPathways = {gly, tca, PPP, ex, pyr, purine, pyrimidine, nuc, glut, ureaCycle, proline, lipids};
                namesPathways = ["Glycolysis/gluconeogenesis", "Citric acid cycle", "Pentose phosphate pathway", "Exchange/demand reaction", "Pyruvate metabolism",...
                                  "Purine metabolism", "Pyrimidine metabolism", "Nucleotide interconversion", "Glutamate metabolism", "Urea cycle",...
                                  "Arginine and proline metabolism", "Lipid metabolism"];
end
