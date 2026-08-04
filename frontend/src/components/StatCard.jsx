import React from 'react';
import Icon from './Icon';

// const StatCard = ({ title, value, subtext, icon, color = "text-brand-600" }) => (
//     <div className="bg-white p-5 rounded-xl border border-slate-200 shadow-sm hover:shadow-md hover:border-brand-300 transition-all duration-200 flex flex-col justify-between h-full">
//         <div className="flex justify-between items-start">
//             <h3 className="text-xs font-bold text-slate-500 uppercase tracking-wider">{title}</h3>
//             {icon && <div className="p-1.5 bg-slate-50 rounded text-slate-400"><Icon name={icon} size={16} /></div>}
//         </div>
//         <div className="mt-3">
//             <div className={`text-2xl font-bold ${color}`}>{value}</div>
//             {subtext && <p className="mt-1 text-xs text-slate-400 font-medium">{subtext}</p>}
//         </div>
//     </div>
// );
const StatCard = ({ title, value, subtext, icon, color = "text-brand-600", onClick, tooltip }) => (                                                      
        <div                                                                                                                                                 
            onClick={onClick}                                                                                                                                
            title={tooltip}                                                                                                                                  
            // Include cursor-pointer and hover:bg-slate-50 only if onClick is provided                                                                              
            className={`bg-white p-5 rounded-xl border border-slate-200 shadow-sm transition-all duration-200 flex flex-col justify-between h-full           
  hover:shadow-md hover:border-brand-300 ${onClick ? 'cursor-pointer hover:bg-slate-50' : ''}`}                                                              
        >                                                                                                                                                    
            <div className="flex justify-between items-start">                                                                                               
                <h3 className="text-xs font-bold text-slate-500 uppercase tracking-wider">{title}</h3>                                                       
                {icon && <div className="p-1.5 bg-slate-50 rounded text-slate-400"><Icon name={icon} size={16} /></div>}                                     
            </div>                                                                                                                                           
            <div className="mt-3">                                                                                                                           
                <div className={`text-2xl font-bold ${color}`}>{value}</div>                                                                                 
                {subtext && <p className="mt-1 text-xs text-slate-400 font-medium">{subtext}</p>}                                                            
            </div>                                                                                                                                           
        </div>                                                                                                                                               
    );

export default StatCard;
