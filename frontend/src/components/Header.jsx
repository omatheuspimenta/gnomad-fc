// frontend/src/components/Header.jsx
import React from 'react';
import Icon from './Icon';
import IconLogo from './IconLogo'; // Import your new component

const Header = ({ activeTab, setActiveTab }) => {
    return (
        <header className="bg-[#0f172a] text-white border-b border-slate-800 sticky top-0 z-50 shadow-lg">
            <div className="max-w-7xl mx-auto px-4 pt-4">
                <div className="flex items-center justify-between pb-4">
                    <div className="flex items-center gap-4">

                        {/* <IconLogo className="h-12 w-auto text-white" /> */}
                        {/* <IconLogo 
                            className="h-12 w-auto text-white cursor-pointer" 
                            onClick={() => setActiveTab('home')}
                            title="Home"
                        /> */}
                        <div 
                            onClick={() => setActiveTab('home')}
                            title="Home"
                            className="cursor-pointer flex items-center"
                                    >
                            <IconLogo className="h-12 w-auto text-white" />
                        </div>
                        
                        <div className="h-8 w-px bg-slate-700 mx-1 hidden sm:block"></div>

                        <div>
                            <h1 className="text-xl font-bold tracking-tight leading-none">aPRoVAR</h1>
                            <p className="text-xs text-slate-400 mt-1 font-light tracking-wide uppercase">Arquivo Paranaense Online de Variantes Genéticas</p>
                        </div>
                    </div>
                    <div className="hidden md:flex items-center gap-4 text-sm text-slate-400">
                        <span className="flex items-center gap-1"><Icon name="database" size={14} /> Connected</span>
                    </div>
                </div>

                {/* Navigation Tabs */}
                <div className="flex items-center gap-8 text-sm font-medium">
                    <button 
                        onClick={() => setActiveTab('home')}
                        className={`pb-3 border-b-2 transition-colors ${activeTab === 'home' ? 'border-brand-500 text-white' : 'border-transparent text-slate-400 hover:text-white'}`}
                    >
                        Home
                    </button>
                    <button 
                        onClick={() => setActiveTab('about')}
                        className={`pb-3 border-b-2 transition-colors ${activeTab === 'about' ? 'border-brand-500 text-white' : 'border-transparent text-slate-400 hover:text-white'}`}
                    >
                        About
                    </button>
                    <button 
                        onClick={() => setActiveTab('terms')}
                        className={`pb-3 border-b-2 transition-colors ${activeTab === 'terms' ? 'border-brand-500 text-white' : 'border-transparent text-slate-400 hover:text-white'}`}
                    >
                        Terms
                    </button>
                    <button 
                        onClick={() => setActiveTab('download')}
                        className={`pb-3 border-b-2 transition-colors ${activeTab === 'download' ? 'border-brand-500 text-white' : 'border-transparent text-slate-400 hover:text-white'}`}
                    >
                        Download
                    </button>
                </div>
            </div>
        </header>
    );
};

export default Header;