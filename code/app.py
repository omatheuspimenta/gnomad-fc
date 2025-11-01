import streamlit as st
import hail as hl
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from gnomad_toolbox.filtering.vep import get_gene_intervals
import re

# Configure Streamlit page
st.set_page_config(page_title="gnomAD Variant Browser", layout="wide")

# helper function to extract info field
def first_or_default(expr, default):
    return hl.if_else(hl.len(expr) > 0, expr[0], default)

@st.cache_resource
def init_hail():
    """Initialize Hail once and cache it"""
    hl.init(quiet=True, log='../logs/hail')
    return True

@st.cache_data
def load_gene_list():
    """Load gene symbols for autocomplete"""
    try:
        # Get gene intervals from gnomad-toolbox
        genes = get_gene_intervals(reference='GRCh38')
        gene_df = genes.to_pandas()
        return sorted(gene_df['gene_name'].unique().tolist())
    except:
        return []

@st.cache_data
def load_variant_data(mt_path):
    """Load and process matrix table"""
    mt = hl.read_matrix_table(mt_path)
    variants = mt.rows()
    
    # Select relevant fields from the info struct
    variants = variants.select(
        # Locus information
        chrom=variants.locus.contig,
        pos=variants.locus.position,
        ref=variants.alleles[0],
        alt=hl.delimit(variants.alleles[1:], ','),  # Join alt alleles if multiple
        rsid=variants.rsid,

        # Basic allele statistics (from gnomAD info field)
        AC=first_or_default(variants.info.AC, 0),
        AN=variants.info.AN,
        AF=first_or_default(variants.info.AF, 0.0),

        # Genetic ancestry group with max AF
        grpmax=first_or_default(variants.info.grpmax, ''),

        # Filtering Allele Frequency (95% confidence)
        faf95_max=first_or_default(variants.info.fafmax_faf95_max, 0.0),
        faf95_max_gen_anc=first_or_default(variants.info.fafmax_faf95_max_gen_anc, ''),

        # Number of homozygotes
        nhomalt=first_or_default(variants.info.nhomalt, 0),

        # Sex-specific counts
        AC_XX=first_or_default(variants.info.AC_XX, 0),
        AN_XX=variants.info.AN_XX,
        AF_XX=first_or_default(variants.info.AF_XX, 0.0),
        nhomalt_XX=first_or_default(variants.info.nhomalt_XX, 0),

        AC_XY=first_or_default(variants.info.AC_XY, 0),
        AN_XY=variants.info.AN_XY,
        AF_XY=first_or_default(variants.info.AF_XY, 0.0),
        nhomalt_XY=first_or_default(variants.info.nhomalt_XY, 0),

        # Population-specific (examples)
        AC_afr=first_or_default(variants.info.AC_afr, 0),
        AF_afr=first_or_default(variants.info.AF_afr, 0.0),
        AN_afr=variants.info.AN_afr,

        AC_amr=first_or_default(variants.info.AC_amr, 0),
        AF_amr=first_or_default(variants.info.AF_amr, 0.0),
        AN_amr=variants.info.AN_amr,

        AC_eas=first_or_default(variants.info.AC_eas, 0),
        AF_eas=first_or_default(variants.info.AF_eas, 0.0),
        AN_eas=variants.info.AN_eas,

        AC_nfe=first_or_default(variants.info.AC_nfe, 0),
        AF_nfe=first_or_default(variants.info.AF_nfe, 0.0),
        AN_nfe=variants.info.AN_nfe,

        AC_sas=first_or_default(variants.info.AC_sas, 0),
        AF_sas=first_or_default(variants.info.AF_sas, 0.0),
        AN_sas=variants.info.AN_sas,

        # Prediction scores
        cadd_phred=variants.info.cadd_phred,
        revel_max=variants.info.revel_max,
        spliceai_ds_max=variants.info.spliceai_ds_max,

        # Quality metrics
        qual=variants.qual,
        filters=hl.delimit(variants.filters, ','),
    )

    
    return variants

def filter_by_gene(variants_ht, gene_name):
    """Filter variants by gene name using VEP annotations"""
    # Filter using VEP transcript consequences
    filtered = variants_ht.filter(
        hl.any(lambda x: x.contains(f"|{gene_name}|"), variants_ht.vep)
    )
    return filtered

def parse_variant_string(variant_str):
    """Parse variant string like chr1-12345-A-T or rs123456"""
    if variant_str.startswith('rs'):
        return None, None, None, None, variant_str
    
    # Try chr-pos-ref-alt format
    match = re.match(r'chr?(\w+)[:-](\d+)[:-]([ACGT]+)[:-]([ACGT]+)', variant_str, re.IGNORECASE)
    if match:
        chrom, pos, ref, alt = match.groups()
        return f"chr{chrom}", int(pos), ref.upper(), alt.upper(), None
    
    return None, None, None, None, None

# Initialize
init_hail()

# Sidebar
with st.sidebar:
    st.header("⚙️ Configuration")
    mt_path = st.text_input("Matrix Table Path", value="../data/merged.mt")
    
    if st.button("🔄 Load/Reload Data"):
        with st.spinner("Loading variant data..."):
            try:
                st.session_state.variants_ht = load_variant_data(mt_path)
                st.session_state.gene_list = load_gene_list()
                st.success("✅ Data loaded!")
            except Exception as e:
                st.error(f"Error loading data: {e}")
    
    st.markdown("---")
    st.markdown("### 📖 Help")
    st.markdown("""
    **Search by:**
    - **Gene**: Enter gene symbol (e.g., BRCA1)
    - **Variant**: Enter variant ID
        - rsID: rs123456
        - Position: chr1-12345-A-T
    """)

# Main content
st.title("🧬 gnomAD Variant Browser")
st.markdown("Interactive browser for gnomAD exome variant data")

if 'variants_ht' not in st.session_state:
    st.info("👈 Please load data using the sidebar")
    st.stop()

variants_ht = st.session_state.variants_ht

# Search section
st.markdown("---")
col1, col2, col3 = st.columns([2, 2, 1])

with col1:
    search_type = st.radio("Search by:", ["Gene", "Variant", "Region"], horizontal=True)

with col2:
    if search_type == "Gene":
        gene_list = st.session_state.get('gene_list', [])
        if gene_list:
            search_query = st.selectbox("Select Gene:", [""] + gene_list)
        else:
            search_query = st.text_input("Enter Gene Symbol:", placeholder="e.g., BRCA1")
    elif search_type == "Variant":
        search_query = st.text_input("Enter Variant:", placeholder="e.g., rs123456 or chr1-12345-A-T")
    else:
        search_query = st.text_input("Enter Region:", placeholder="e.g., chr1:1000000-2000000")

with col3:
    st.write("")  # Spacing
    st.write("")  # Spacing
    search_button = st.button("🔍 Search", type="primary", use_container_width=True)

# Execute search
filtered_ht = None
if search_button and search_query:
    with st.spinner("Searching..."):
        try:
            if search_type == "Gene":
                filtered_ht = filter_by_gene(variants_ht, search_query.upper())
                search_desc = f"Gene: {search_query}"
            
            elif search_type == "Variant":
                chrom, pos, ref, alt, rsid = parse_variant_string(search_query)
                if rsid:
                    filtered_ht = variants_ht.filter(variants_ht.rsid == rsid)
                    search_desc = f"rsID: {rsid}"
                elif chrom and pos:
                    filtered_ht = variants_ht.filter(
                        (variants_ht.chrom == chrom) &
                        (variants_ht.pos == pos)
                    )
                    if ref and alt:
                        filtered_ht = filtered_ht.filter(
                            (filtered_ht.ref == ref) &
                            (filtered_ht.alt == alt)
                        )
                    search_desc = f"Variant: {search_query}"
                else:
                    st.error("Invalid variant format")
                    st.stop()
            
            else:  # Region
                match = re.match(r'chr?(\w+)[:-](\d+)-(\d+)', search_query, re.IGNORECASE)
                if match:
                    chrom, start, end = match.groups()
                    chrom = f"chr{chrom}"
                    filtered_ht = variants_ht.filter(
                        (variants_ht.chrom == chrom) &
                        (variants_ht.pos >= int(start)) &
                        (variants_ht.pos <= int(end))
                    )
                    search_desc = f"Region: {search_query}"
                else:
                    st.error("Invalid region format")
                    st.stop()
            
            # Convert to pandas
            df = filtered_ht.to_pandas()
            
            if len(df) == 0:
                st.warning(f"No variants found for {search_desc}")
            else:
                st.success(f"Found {len(df):,} variants for {search_desc}")
                st.session_state.current_df = df
                st.session_state.search_desc = search_desc
        
        except Exception as e:
            st.error(f"Search error: {e}")

# Display results
if 'current_df' in st.session_state:
    df = st.session_state.current_df
    
    # Summary metrics
    st.markdown("---")
    st.subheader(f"📊 Summary: {st.session_state.search_desc}")
    
    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.metric("Variants", f"{len(df):,}")
    with col2:
        st.metric("Mean AC", f"{df['AC'].mean():.2f}")
    with col3:
        st.metric("Mean AF", f"{df['AF'].mean():.6f}")
    with col4:
        st.metric("Mean FAF95", f"{df['faf95_max'].mean():.6f}")
    with col5:
        st.metric("Total Homozygotes", f"{df['nhomalt'].sum():,}")
    
    # Filters
    st.markdown("---")
    with st.expander("🔧 Apply Filters", expanded=False):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            af_min = st.number_input("Min AF", min_value=0.0, max_value=1.0, value=0.0, format="%.6f")
            af_max = st.number_input("Max AF", min_value=0.0, max_value=1.0, value=1.0, format="%.6f")
        
        with col2:
            ac_min = st.number_input("Min AC", min_value=0, value=0)
            ac_max = st.number_input("Max AC", min_value=0, value=int(df['AC'].max()))
        
        with col3:
            min_cadd = st.number_input("Min CADD", min_value=0.0, value=0.0)
            min_revel = st.number_input("Min REVEL", min_value=0.0, max_value=1.0, value=0.0)
        
        # Apply filters
        mask = (
            (df['AF'] >= af_min) & (df['AF'] <= af_max) &
            (df['AC'] >= ac_min) & (df['AC'] <= ac_max)
        )
        if min_cadd > 0:
            mask &= (df['cadd_phred'] >= min_cadd)
        if min_revel > 0:
            mask &= (df['revel_max'] >= min_revel)
        
        df = df[mask]
        st.info(f"After filtering: {len(df):,} variants")
    
    # Visualizations
    st.markdown("---")
    tab1, tab2, tab3, tab4 = st.tabs(["📈 AF Distribution", "🌍 Population AF", "🔬 Predictions", "📋 Variant Table"])
    
    with tab1:
        col1, col2 = st.columns(2)
        with col1:
            fig = px.histogram(df, x='AF', nbins=50, title="Allele Frequency Distribution")
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            fig = px.scatter(df, x='AC', y='AF', hover_data=['rsid', 'chrom', 'pos'],
                           title="Allele Count vs Frequency", opacity=0.6)
            st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        pop_cols = ['AF_afr', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_sas']
        pop_names = ['African', 'Latino', 'East Asian', 'European', 'South Asian']
        
        # Average AF by population
        pop_means = [df[col].mean() for col in pop_cols]
        fig = go.Figure(data=[go.Bar(x=pop_names, y=pop_means)])
        fig.update_layout(title="Mean Allele Frequency by Population", 
                         xaxis_title="Population", yaxis_title="Mean AF")
        st.plotly_chart(fig, use_container_width=True)
        
        # Heatmap for variants
        if len(df) <= 100:
            pop_data = df[pop_cols].fillna(0)
            fig = px.imshow(pop_data.T, labels=dict(x="Variant Index", y="Population", color="AF"),
                          y=pop_names, aspect="auto", title="Population-specific Allele Frequencies")
            st.plotly_chart(fig, use_container_width=True)
    
    with tab3:
        col1, col2 = st.columns(2)
        with col1:
            fig = px.scatter(df, x='cadd_phred', y='revel_max', 
                           hover_data=['rsid', 'chrom', 'pos', 'AF'],
                           title="CADD vs REVEL Scores", opacity=0.6)
            st.plotly_chart(fig, use_container_width=True)
        with col2:
            fig = px.histogram(df, x='spliceai_ds_max', nbins=30,
                             title="SpliceAI Score Distribution")
            st.plotly_chart(fig, use_container_width=True)
    
    with tab4:
        # Display columns
        display_cols = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'AC', 'AN', 'AF', 
                       'faf95_max', 'nhomalt', 'grpmax', 'cadd_phred', 'revel_max']
        
        st.dataframe(
            df[display_cols].style.format({
                'AF': '{:.6f}',
                'faf95_max': '{:.6f}',
                'cadd_phred': '{:.2f}',
                'revel_max': '{:.3f}'
            }),
            use_container_width=True,
            height=400
        )
        
        # Download
        csv = df.to_csv(index=False)
        st.download_button(
            "📥 Download Results (CSV)",
            data=csv,
            file_name=f"gnomad_variants_{st.session_state.search_desc.replace(' ', '_')}.csv",
            mime="text/csv"
        )