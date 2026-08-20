import React, { useState, useEffect } from 'react';
import Header from './components/Header';
import SearchHero from './components/SearchHero';
import StatCard from './components/StatCard';
import Icon from './components/Icon';
import VariantTable from './components/VariantTable';
import FrequencyPie from './components/Charts/FrequencyPie';
import PopulationBar from './components/Charts/PopulationBar';
import VariantTypeBar from './components/Charts/VariantTypeBar';
import QualityHistogram from './components/Charts/QualityHistogram';
import ConservationScatter from './components/Charts/ConservationScatter';
import ConservationBar from './components/Charts/ConservationBar';
import AboutThisProject from './components/AboutThisProject';
import { useSearch } from './hooks/useSearch';

function App() {
    const {
        rawData, filteredVariants, loading, error,
        filters, setFilters, handleSearch, stats,
        currentPage, totalPages, changePage, totalVariants,
        currentSearchType, currentQuery
    } = useSearch();

    const [activeTab, setActiveTab] = useState('frequency');
    const [showFilters, setShowFilters] = useState(false);
    const [mainTab, setMainTab] = useState('home');


    useEffect(() => {
            // The condition below ensures that the alert will only appear if the user has already performed a search (i.e., if there is data to be lost).
            // If you want the alert to appear ALWAYS (even with an empty screen), just remove the "if" below.
            if (!rawData || rawData.length === 0) {
                return;
            }
  
            const handleBeforeUnload = (event) => {
                // O preventDefault and the event.returnValue are required for most modern browsers to trigger the confirmation dialog.
                event.preventDefault();
                event.returnValue = '';
                return '';
            };
  
            // Add the "spy" to listen for the beforeunload event, which is triggered when the user tries to close the tab or navigate away.
            window.addEventListener('beforeunload', handleBeforeUnload);
  
            // Remove the "spy" when the component is unmounted to avoid memory leaks
            return () => {
                window.removeEventListener('beforeunload', handleBeforeUnload);
            };
        }, [rawData]); // The array [rawData] makes the React update this rule whenever the search data changes.

    const requestBody = `Hello, I would like to request data for the following search parameters:

Search Type: ${currentSearchType || 'N/A'}
Search Query: ${currentQuery || 'N/A'}

Filters applied:
${JSON.stringify(filters, null, 2)}`;

    return (
        <div className="min-h-screen flex flex-col bg-slate-50 text-slate-900 font-sans">
            <Header activeTab={mainTab} setActiveTab={setMainTab} />

            {mainTab === 'home' && (
                <>
                    <SearchHero onSearch={handleSearch} loading={loading} />

                    <main className="flex-1 max-w-7xl w-full mx-auto px-4 py-8">
                        {error && (
                            <div className="bg-red-50 border border-red-200 text-red-800 p-4 rounded-lg flex items-center gap-3 mb-8 animate-fade-in shadow-sm">
                                <Icon name="alert-triangle" />
                                <p className="font-medium">{error}</p>
                            </div>
                        )}

                        {!rawData && !error && !loading && (
                            <AboutThisProject />
                        )}

                        {rawData && (
                    <div className="space-y-6 animate-fade-in">
                        {stats && (
                            <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-6 gap-4">
                                {/* <StatCard title="Total Variants" value={stats.count.toLocaleString()} icon="layers" /> */}
                                <StatCard                                                                                                                                                
                                    title="Total Variants"                                                                                                                               
                                    value={stats.count.toLocaleString()}                                                                                                                 
                                    icon="layers"                                                                                                                                        
                                    onClick={() => setActiveTab('table')}                                                                                                                
                                    tooltip="Click to view the data table"                                                                                
                                /> 
                                {/* <StatCard title="Variant Types" value={stats.uniqueTypes} icon="git-branch" /> */}
                                <StatCard                                                                                                                                                
                                    title="Variant Types"                                                                                                                                
                                    value={stats.uniqueTypes}                                                                                                                            
                                    icon="git-branch"                                                                                                                                    
                                    onClick={() => setActiveTab('variantTypes')}                                                                                                         
                                    tooltip="Click to view charts and distribution by variant types"                                                                           
                                />
                                {/* <StatCard title="Mean AF" value={stats.meanAF.toExponential(1)} icon="activity" subtext="Global gnomAD" /> */}
                                {/* <StatCard title="Mean MAF aPRoVAR" value={stats.localMeanAF.toFixed(4)} icon="activity" subtext="Mean MAF aPRoVAR" /> */}
                                <StatCard                                                                                                                                                
                                    title="Mean MAF aPRoVAR"                                                                                                                             
                                    value={stats.localMeanAF.toFixed(4)}                                                                                                                 
                                    icon="activity"                                                                                                                                      
                                    subtext="Mean MAF aPRoVAR"                                                                                                                           
                                    tooltip="Distribution by allele frequency - aPRoVAR"
                                />
                                {/* <StatCard title="Mean MAF gnomAD" value={stats.gnomadMeanAF.toFixed(3)} icon="trending-up" /> */}
                                <StatCard                                                                                                                                                
                                    title="Mean MAF gnomAD"                                                                                                                              
                                    value={stats.gnomadMeanAF.toFixed(3)}                                                                                                                
                                    icon="trending-up"                                                                                                                                   
                                    tooltip="Distribution by allele frequency - gnomAD"
                                />
                                {/* <StatCard title="ClinVar" value={stats.clinvarCount} icon="clipboard-check" color="text-purple-600" /> */}
                                <StatCard                                                                                                                                                
                                    title="ClinVar"                                                                                                                                      
                                    value={stats.clinvarCount}                                                                                                                           
                                    icon="clipboard-check"                                                                                                                               
                                    color="text-purple-600" 
                                    tooltip="Variants with clinical significance reported in ClinVar"
                                />
                                {/* <StatCard title="Coverage" value={stats.coverage} icon="pie-chart" color="text-slate-600" /> */}
                                <StatCard 
                                    title="Coverage" 
                                    value={stats.coverage} 
                                    icon="pie-chart" 
                                    color="text-slate-600" 
                                    tooltip="This metric represents the proportion of variants with coverage data."
                                />
                            </div>
                        )}

                        <div className="bg-white rounded-xl border border-slate-200 shadow-sm overflow-hidden">
                            <div className="p-4 border-b border-slate-100 bg-slate-50/50 flex flex-col md:flex-row md:items-center justify-between gap-4">
                                <div className="flex gap-1 bg-slate-200/60 p-1 rounded-lg overflow-x-auto no-scrollbar">
                                    {[
                                        { id: 'frequency', label: 'Frequency', icon: 'bar-chart-2' },
                                        { id: 'populations', label: 'Populations', icon: 'globe' },
                                        { id: 'conservation', label: 'Conservation', icon: 'shield' },
                                        { id: 'variantTypes', label: 'Types', icon: 'dna' },
                                        { id: 'quality', label: 'Quality', icon: 'check-circle' },
                                        { id: 'table', label: 'Data Table', icon: 'table' }
                                    ].map(tab => (
                                        <button
                                            key={tab.id}
                                            onClick={() => setActiveTab(tab.id)}
                                            className={`px-3 py-1.5 text-sm font-medium rounded-md flex items-center gap-2 transition-all whitespace-nowrap ${activeTab === tab.id ? 'bg-white text-brand-700 shadow-sm' : 'text-slate-600 hover:text-slate-900 hover:bg-slate-200/50'}`}
                                        >
                                            <Icon name={tab.icon} size={14} /> {tab.label}
                                        </button>
                                    ))}
                                </div>
                                <div className="flex items-center gap-2">
                                    <button
                                        onClick={() => setShowFilters(!showFilters)}
                                        className={`px-4 py-2 text-sm font-medium rounded-lg border flex items-center gap-2 transition-colors ${showFilters ? 'bg-brand-50 border-brand-200 text-brand-700' : 'bg-white border-slate-200 text-slate-600 hover:bg-slate-50'}`}
                                    >
                                        <Icon name="filter" size={16} /> Filters
                                    </button>
                                    <a 
                                        href={`mailto:tiago.gomes@fiocruz.br?subject=${encodeURIComponent("[DATA REQUEST] gnomAD-tb Search")}&body=${encodeURIComponent(requestBody)}`}
                                        className={`px-4 py-2 text-sm font-medium text-white rounded-lg flex items-center gap-2 shadow-sm ${!filteredVariants.length ? 'bg-slate-300 cursor-not-allowed pointer-events-none' : 'bg-emerald-600 hover:bg-emerald-700'}`}
                                    >
                                        <Icon name="mail" size={16} /> Request
                                    </a>
                                </div>
                            </div>

                            {showFilters && (
                                <div className="p-6 bg-slate-50 border-b border-slate-200 grid grid-cols-1 md:grid-cols-3 gap-6 animate-fade-in">
                                    {/* <div className="space-y-2">
                                        <label className="text-xs font-bold text-slate-500 uppercase">Frequency (AF) Range</label>
                                        <div className="flex items-center gap-2">
                                            <input type="range" min="0" max="1" step="0.0001" value={filters.minAF} onChange={e => setFilters({ ...filters, minAF: parseFloat(e.target.value) })} className="flex-1" />
                                            <span className="text-xs font-mono bg-white border px-2 py-1 rounded">{filters.minAF.toFixed(4)}</span>
                                        </div>
                                    </div>
                                    <div className="space-y-2">
                                        <label className="text-xs font-bold text-slate-500 uppercase">Max Missing (F_MISSING)</label>
                                        <div className="flex items-center gap-2">
                                            <input type="range" min="0" max="1" step="0.01" value={filters.maxFMissing} onChange={e => setFilters({ ...filters, maxFMissing: parseFloat(e.target.value) })} className="flex-1" />
                                            <span className="text-xs font-mono bg-white border px-2 py-1 rounded">{filters.maxFMissing.toFixed(2)}</span>
                                        </div>
                                    </div> */}
                                    <div className="space-y-2">                                                                                                                              
                                        <label className="text-xs font-bold text-slate-500 uppercase">Frequency (AF) Range</label>                                                           
                                        <div className="flex items-center gap-2">                                                                                                            
                                            <input type="range" min="0" max="1" step="0.0001" value={filters.minAF} onChange={e => setFilters({ ...filters, minAF: parseFloat(e.target.value)
                                })} className="flex-1" />                                                                                                                                  
                                            {/* new input tag below */}                                                                                                         
                                            <input                                                                                                                                           
                                                type="number"                                                                                                                                
                                                min="0"                                                                                                                                      
                                                max="1"                                                                                                                                      
                                                step="0.0001"                                                                                                                                
                                                value={filters.minAF}                                                                                                                        
                                                onChange={e => {                                                                                                                             
                                                    // prevents errors in case of user deletes everything in the input (empty field)                                                                       
                                                    const val = e.target.value === '' ? 0 : parseFloat(e.target.value);                                                                      
                                                    setFilters({ ...filters, minAF: val });                                                                                                  
                                                }}                                                                                                                                           
                                                className="text-xs font-mono bg-white border border-slate-300 px-2 py-1 rounded w-24 outline-none focus:ring-2 focus:ring-brand-500          
                                focus:border-brand-500"                                                                                                                                    
                                            />                                                                                                                                               
                                        </div>                                                                                                                                               
                                    </div> 
                                    <div className="space-y-2">
                                        <label className="text-xs font-bold text-slate-500 uppercase">Max Missing (F_MISSING)</label>
                                        <div className="flex items-center gap-2">
                                            <input type="range" min="0" max="1" step="0.01" value={filters.maxFMissing} onChange={e => setFilters({ ...filters, maxFMissing: parseFloat(e.   
                                target.value) })} className="flex-1" />
                                            {/* new input tag include below */}
                                            <input 
                                                type="number" 
                                                min="0" 
                                                max="1" 
                                                step="0.01" 
                                                value={filters.maxFMissing} 
                                                onChange={e => {
                                                    // Prevents errors in case the user deletes everything in the input (empty field)
                                                    const val = e.target.value === '' ? 0 : parseFloat(e.target.value);
                                                    setFilters({ ...filters, maxFMissing: val });
                                                }}
                                                className="text-xs font-mono bg-white border border-slate-300 px-2 py-1 rounded w-24 outline-none focus:ring-2 focus:ring-brand-500          
                                focus:border-brand-500" 
                                            />
                                        </div>
                                    </div>
                                    {['consequence', 'variantType', 'clinvar'].map(f => (
                                        <div key={f} className="space-y-2">
                                            <label className="text-xs font-bold text-slate-500 uppercase">{f.replace(/([A-Z])/g, ' $1').trim()}</label>
                                            <input
                                                type="text"
                                                placeholder="Filter..."
                                                value={filters[f]}
                                                onChange={e => setFilters({ ...filters, [f]: e.target.value })}
                                                className="w-full px-3 py-2 text-sm border border-slate-300 rounded-lg focus:ring-2 focus:ring-brand-500 outline-none"
                                            />
                                        </div>
                                    ))}
                                </div>
                            )}

                            <div className="p-6 min-h-[500px]">
                                {filteredVariants.length === 0 ? (
                                    <div className="text-center py-24 animate-fade-in">
                                        <div className="inline-block p-4 bg-slate-100 rounded-full mb-4 text-slate-400">
                                            <Icon name="search-x" size={48} />
                                        </div>
                                        <h2 className="text-xl font-semibold text-slate-700">No variants found</h2>
                                        <p className="text-slate-500 mt-2">Try adjusting filters or your search query.</p>
                                    </div>
                                ) : (
                                    stats && (
                                        <>
                                            {activeTab === 'frequency' && (
                                                <div className="grid grid-cols-1 lg:grid-cols-3 gap-8 animate-fade-in h-[450px]">
                                                    <div className="lg:col-span-2 bg-white border border-slate-100 rounded-xl p-4 shadow-sm">
                                                        <h4 className="text-sm font-bold text-slate-700 mb-4">Allele Frequency Distribution</h4>
                                                        <FrequencyPie data={stats.pieData} />
                                                    </div>
                                                    <div className="bg-brand-50/50 border border-brand-100 rounded-xl p-6">
                                                        <h4 className="text-sm font-bold text-brand-900 uppercase mb-4">Quick Insights</h4>
                                                        <ul className="space-y-4 text-sm text-brand-800">
                                                            <li className="flex gap-3">
                                                                <Icon name="info" size={18} className="text-brand-500 shrink-0" />
                                                                <span>The majority of variants ({stats.localPieData?.[0]?.value}) fall into the <strong>{stats.localPieData?.[0]?.name}</strong> category in Paraná.</span>
                                                            </li>
                                                            <li className="flex gap-3">
                                                                <Icon name="alert-circle" size={18} className="text-purple-500 shrink-0" />
                                                                <span><strong>{stats.clinvarCount}</strong> variants have clinical significance reported in ClinVar.</span>
                                                            </li>
                                                            <li className="flex gap-3">
                                                                <Icon name="trending-up" size={18} className="text-emerald-500 shrink-0" />
                                                                <span>Mean Allele Frequency observed globally is <strong>{(stats.gnomadMeanAF * 100).toFixed(2)}%</strong>.</span>
                                                            </li>
                                                        </ul>
                                                    </div>
                                                </div>
                                            )}

                                            {activeTab === 'populations' && (
                                                <div className="h-[450px] bg-white border border-slate-100 rounded-xl p-6 shadow-sm animate-fade-in">
                                                    <h4 className="text-sm font-bold text-slate-700 mb-6">Global Population Frequencies (Mean AF)</h4>
                                                    <PopulationBar data={stats.popData} />
                                                </div>
                                            )}

                                            {activeTab === 'variantTypes' && (
                                                <div className="grid grid-cols-1 lg:grid-cols-3 gap-8 animate-fade-in h-[450px]">
                                                    <div className="lg:col-span-2 bg-white border border-slate-100 rounded-xl p-6 shadow-sm">
                                                        <h4 className="text-sm font-bold text-slate-700 mb-2">Variant Type Distribution</h4>
                                                        <VariantTypeBar data={stats.variantTypeData} />
                                                    </div>
                                                    <div className="bg-emerald-50/50 border border-emerald-100 rounded-xl p-6 overflow-y-auto custom-scroll">
                                                        <h4 className="text-sm font-bold text-emerald-900 uppercase mb-4">Breakdown</h4>
                                                        <div className="space-y-3">
                                                            {stats.variantTypeData.map((type, i) => (
                                                                <div key={i} className="flex justify-between items-center bg-white p-3 rounded-lg shadow-sm">
                                                                    <span className="text-sm text-slate-600 font-medium">{type.name}</span>
                                                                    <span className="text-sm font-bold text-emerald-600">{type.value}</span>
                                                                </div>
                                                            ))}
                                                        </div>
                                                    </div>
                                                </div>
                                            )}

                                            {activeTab === 'quality' && (
                                                <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 animate-fade-in h-[450px]">
                                                    <div className="bg-white border border-slate-100 rounded-xl p-6 shadow-sm">
                                                        <h4 className="text-sm font-bold text-slate-700 mb-4">Quality Score Histogram</h4>
                                                        <QualityHistogram data={stats.qualityDist} />
                                                    </div>
                                                    <div className="bg-amber-50/50 border border-amber-100 rounded-xl p-6">
                                                        <h4 className="text-sm font-bold text-amber-900 uppercase mb-4">Quality Analysis</h4>
                                                        <p className="text-sm text-slate-600 mb-4">
                                                            Quality scores indicate the confidence in the variant call. Higher scores (right side) generally indicate higher confidence and better sequencing depth.
                                                        </p>
                                                        <div className="space-y-2">
                                                            {stats.qualityDist.map((d, i) => (
                                                                <div key={i} className="flex items-center gap-2">
                                                                    <div className="w-24 text-xs font-bold text-slate-500">{d.name}</div>
                                                                    <div className="flex-1 h-2 bg-slate-200 rounded-full overflow-hidden">
                                                                        <div className="h-full bg-amber-400" style={{ width: `${Math.max(...stats.qualityDist.map(x => x.value)) > 0 ? (d.value / Math.max(...stats.qualityDist.map(x => x.value))) * 100 : 0}%` }}></div>
                                                                    </div>
                                                                    <div className="w-12 text-xs text-right font-mono">{d.value}</div>
                                                                </div>
                                                            ))}
                                                        </div>
                                                    </div>
                                                </div>
                                            )}

                                            {activeTab === 'conservation' && (
                                                <div className="grid grid-cols-1 lg:grid-cols-2 gap-8 animate-fade-in h-[450px]">
                                                    <div className="bg-white border border-slate-100 rounded-xl p-6 shadow-sm">
                                                        <h4 className="text-sm font-bold text-slate-700 mb-2">PhyloP vs Allele Frequency</h4>
                                                        <ConservationScatter data={stats.scatterData} />
                                                    </div>
                                                    <div className="bg-white border border-slate-100 rounded-xl p-6 shadow-sm">
                                                        <h4 className="text-sm font-bold text-slate-700 mb-2">Conservation Score Averages</h4>
                                                        <ConservationBar data={stats.conservationData} />
                                                    </div>
                                                </div>
                                            )}

                                            {activeTab === 'table' && (
                                                <div className="animate-fade-in">
                                                    <VariantTable
                                                        variants={filteredVariants}
                                                        onVariantClick={(vid) => handleSearch('variant', vid)}
                                                        currentPage={currentPage}
                                                        totalPages={totalPages}
                                                        onPageChange={changePage}
                                                        totalVariants={totalVariants}
                                                    />
                                                    <div className="mt-4 text-right text-xs text-slate-400">
                                                        Showing {filteredVariants.length} variants
                                                    </div>
                                                </div>
                                            )}
                                        </>
                                    )
                                )}
                            </div>
                        </div>
                    </div>
                )}
            </main>
                </>
            )}

            {mainTab === 'about' && (
                <main className="flex-1 max-w-7xl w-full mx-auto px-4 py-8 animate-fade-in">
                    {/* <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 max-w-4xl mx-auto block">
                        <h2 className="text-xl font-bold text-slate-800 mb-4">About</h2>
                        <p className="text-slate-600 mb-4">
                            Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
                        </p>
                    </div> */}
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
                </main>
            )}

            {mainTab === 'terms' && (
                <main className="flex-1 max-w-7xl w-full mx-auto px-4 py-8 animate-fade-in">
                    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 max-w-4xl mx-auto block">
                        <h2 className="text-xl font-bold text-slate-800 mb-4">Terms of Use</h2>
                        <p className="text-slate-600 mb-4">
                            Lorem ipsum dolor sit amet, consectetur adipiscing elit. Nullam in dui mauris. 
                            Vivamus hendrerit arcu sed erat molestie vehicula. Sed auctor neque eu tellus rhoncus ut eleifend nibh porttitor.
                            Teste de texto
                        </p>
                    </div>
                </main>
            )}

            {mainTab === 'download' && (
                <main className="flex-1 max-w-7xl w-full mx-auto px-4 py-8 animate-fade-in">
                    <div className="bg-white rounded-xl shadow-sm border border-slate-200 p-8 max-w-4xl mx-auto block">
                        <h2 className="text-xl font-bold text-slate-800 mb-4">Download Data</h2>
                        <p className="text-slate-600 mb-6 font-medium">
                            Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
                        </p>
                        <div className="bg-slate-50 p-6 rounded-lg border border-slate-200">
                            <h3 className="text-lg font-semibold text-slate-800 mb-3">Available Datasets</h3>
                            <ul className="list-disc list-inside text-slate-600 space-y-2 ml-2 mb-6">
                                <li>Item 1</li>
                                <li>Item 2</li>
                                <li>Item 3</li>
                            </ul>
                            <a 
                                href="mailto:tiago.gomes@fiocruz.br?subject=[DATA REQUEST] aPRoVAR"
                                className="inline-block bg-brand-600 hover:bg-brand-700 text-white px-6 py-2 rounded-lg font-medium transition-colors"
                            >
                                Request Data
                            </a>
                        </div>
                    </div>
                </main>
            )}
        </div>
    );
}

export default App;
