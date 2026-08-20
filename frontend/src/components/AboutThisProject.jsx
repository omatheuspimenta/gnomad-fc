import React from 'react';

const AboutThisProject = () => {
    return (
        <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 mt-8 animate-fade-in max-w-4xl mx-auto">
            <h2 className="text-xl font-bold text-slate-800 mb-4">About This Project</h2>
            <p className="text-slate-600 mb-6">
                This is a comprehensive regional database of germline variants from Paraná. This resource provides information to researchers and health professionals for variant interpretation. 
                <br /><br />
                aPRoVAR is based on whole-exome sequencing data of 1,010 individuals from Paraná, Brazil. Allele frequency estimates are adjusted by phenotype group. Variants with allele frequency equal to 0 correspond to those which all alternative genotypes occurred in individuals affected by one of the cohort phenotypes (COVID-19, breast cancer and sepsis). 
            </p>
            
            <div className="grid grid-cols-1 md:grid-cols-2 gap-8 mb-6">
                <div>
                    <h3 className="text-lg font-semibold text-slate-800 mb-3">How to Use</h3>
                    <p className="text-slate-600 mb-3">
                        Use the search bar above to query our database:
                    </p>
                    <ul className="list-disc list-inside text-slate-600 space-y-2 ml-2">
                        <li><strong>By Gene:</strong> Enter a gene symbol (e.g., <i>TP53</i>)</li>
                        <li><strong>By Variant:</strong> Enter a specific variant (e.g., 17:7670699-C-A or rsID)</li>
                        <li><strong>By Region:</strong> Enter a genomic region (e.g., 17:7670000-7671000)</li>
                    </ul>
                </div>
                
                <div>
                    <h3 className="text-lg font-semibold text-slate-800 mb-3">Basic Features</h3>
                    <ul className="list-disc list-inside text-slate-600 space-y-2 ml-2">
                        <li>Search across multiple variant types and genomic regions.</li>
                        <li>View detailed genetic variant annotations.</li>
                        <li>Filter results by functional consequence, variant type, pathogenicity annotation and genotype missingness rate.</li>
                        <li>Request filtered variant data for offline analysis.</li>
                    </ul>
                </div>
            </div>

            <hr className="my-6 border-slate-200" />
            
            <h3 className="text-lg font-semibold text-slate-800 mb-2">The Team</h3>

            <div className="text-slate-600 mb-6">
            <p>
                <span className="font-semibold">Faculty:</span>{" "}
                Prof. Fabio Passetti, PhD; Prof. Helisson Faoro, PhD; and Prof. Hellen
                Geremias Gatica Santos, PhD
            </p>

            <p>
                <span className="font-semibold">Cohort provided by:</span>{" "}
                XXX
            </p>

            <p>
                <span className="font-semibold">aPRoVAR leaders:</span>{" "}
                Marco A. Campanário, MSc; Bruno J. do Nascimento, MSc; and Tiago M. F. F.
                Gomes, PhD
            </p>

            <p>
                <span className="font-semibold">Technical and scientific support:</span>{" "}
                Matheus Henrique Pimenta-Zanon, PhD; Eduardo Martin Tarazona Santos, PhD;
                Michel S. Naslavsky, PhD
            </p>
            </div>

            <hr className="my-6 border-slate-200" />
            
            <h3 className="text-lg font-semibold text-slate-800 mb-2">Citation</h3>
            <p className="text-slate-600 italic mb-6">
            <span className="not-italic">
                By using this resource, you agree to cite our paper:
            </span>
            <br />
            <span className="font-semibold not-italic">
                "aPRoVAR database: a public online resource of 1,010 exomes from an admixed population in Paraná, Brazil"
            </span>
            <span className="not-italic">
                {" "}by Campanário & Janke
            </span>{" "}
            <i>et al.</i>
            <span className="not-italic">
                , 2026 (Citation Details Soon).
            </span>
            </p>

            <div className="mt-8 pt-6 border-t border-slate-100 text-xs text-slate-500">
                <p>
                    <strong>Disclaimer:</strong> The data provided in this project is for academic and research purposes only and may not be used for commercial purposes.
                </p>
            </div>
        </div>
    );
};

export default AboutThisProject;
