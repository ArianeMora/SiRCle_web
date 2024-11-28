import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scircm import SciRCM
from matplotlib import rcParams
from threading import Thread
from sciviso import Sankeyplot
import plotly.express as px

###################################
from streamlit.runtime.scriptrunner import add_script_run_ctx
from functionforDownloadButtons import download_button
###################################

# """
# The base of this app was developed from:
# https://share.streamlit.io/streamlit/example-app-csv-wrangler/
# """


def _max_width_():
    max_width_str = f"max-width: 1800px;"
    st.markdown(
        f"""
    <style>
    .reportview-container .main .block-container{{
        {max_width_str}
    }}
    </style>    
    """,
        unsafe_allow_html=True,
    )

st.title("Signature Regulatory Clustering Model (Beta mode)")

background_method = st.text_input('Background Method', 'P|(M&R)')
c1, c2, c3 = st.columns([4, 4, 4])

df = None
with c1:

    methylation = st.file_uploader(
        "DNA methylation file",
        key="1",
        help="To activate 'wide mode', go to the hamburger menu > Settings > turn on 'wide mode'",
    )
with c2:
    rna = st.file_uploader(
        "DE RNA file",
        key="2",
        help="To activate 'wide mode', go to the hamburger menu > Settings > turn on 'wide mode'",
    )
with c3:
    protein = st.file_uploader(
        "DE protein file",
        key="3",
        help="To activate 'wide mode', go to the hamburger menu > Settings > turn on 'wide mode'",
    )

with c1:
    if methylation is not None:
        methylation = pd.read_csv(methylation)

        meth_logfc_col = st.selectbox('Methylation diff', list(methylation.columns))
        meth_padj_col = st.selectbox('Methylation padj', list(methylation.columns))
        meth_logfc = st.text_input('Methylation diff cutoff', 0.1)
        meth_padj = st.text_input('Methylation padj cutoff', 0.05)
with c2:
    if rna is not None:
        rna = pd.read_csv(rna)
        rna_logfc_col = st.selectbox('RNA diff column', list(rna.columns))
        rna_padj_col = st.selectbox('RNA padj column', list(rna.columns))
        rna_logfc = st.text_input('RNA logFC cutoff', 1.0)
        rna_padj = st.text_input('RNA padj cutoff', 0.05)
with c3:
    if protein is not None:
        protein = pd.read_csv(protein)
        prot_logfc_col = st.selectbox('Protein diff column', list(protein.columns))
        prot_padj_col = st.selectbox('Protein padj column', list(protein.columns))
        prot_logfc = st.text_input('Protein logFC cutoff', 0.5)
        prot_padj_cutoff = st.text_input('Protein padj cutoff', 0.05)

padd3, c0, padd4 = st.columns([1, 6, 1])


def rcm_runner():
    st.info(f'Running SiRCle using your files! Results will appear below shortly.')

    rcm = SciRCM(methylation, rna, protein,
                 rna_logfc_col, rna_padj_col, meth_logfc_col, meth_padj_col, prot_logfc_col, prot_padj_col,
                 gene_id, sep=',',
                 rna_padj_cutoff=float(rna_padj),
                 prot_padj_cutoff=float(prot_padj_cutoff),
                 meth_padj_cutoff=float(meth_padj),
                 rna_logfc_cutoff=float(rna_logfc),
                 prot_logfc_cutoff=float(prot_logfc),
                 meth_diff_cutoff=float(meth_logfc),
                 output_dir='.',
                 non_coding_genes=['None'],
                 output_filename='SiRCLe',
                 bg_type=background_method
             )

    rcm.run()
    df = rcm.get_df()

    # Reduce the size by removing the background data
    plot_df = df[df['RG2_Changes_filtered'] != 'None']
    plot_df = plot_df[plot_df['RG2_Changes_filtered'] != 'Not-Background']
    plot_df = plot_df[plot_df['RG2_Changes_filtered'] != '']
    plot_df = plot_df[plot_df['RG2_Changes_filtered'] != None]

    download_button(
        plot_df,
        "SiRCle.csv",
        "Download SiRCle clustering to CSV",
    )

    sk = Sankeyplot(plot_df)
    fig = sk.plot(columns=['Methylation', 'RNA', 'Protein', 'RG2_Changes',
                           'RG2_Changes_filtered'], colour_col='Protein')
    st.plotly_chart(fig)

    # Define the color map
    colour_map = {
        'MDS': '#d8419b', 'MDS_TMDE': '#e585c0', 'MDS_ncRNA': '#d880b4',
        'MDE': '#6aaf44', 'MDE_TMDS': '#0e8e6d', 'MDE_ncRNA': '#9edb77',
        'TMDE': '#fe2323', 'TMDS': '#2952ff',
        'TPDE': '#e68e25', 'TPDE_TMDS': '#844c0f',
        'TPDS': '#462d76', 'TPDS_TMDE': '#9b29b7'
    }

    # Order of labels
    rcm_labels = ["MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS", "TPDS", "TPDS_TMDE"]
    colours = [colour_map[label] for label in rcm_labels]

    # Create the Plotly histogram
    fig = px.histogram(
        data_frame=plot_df,
        x='RG2_Changes_filtered',
        title='SiRCle group: RG2_Changes_filtered',
        category_orders={'RG2_Changes_filtered': rcm_labels},
        color_discrete_sequence=colours
    )

    # Customize layout
    fig.update_layout(
        xaxis_title=None,
        yaxis_title="Count",
        xaxis_tickangle=45,
        title_x=0.5,
        height=400
    )

    # Streamlit app
    st.title("SiRCle Groups Visualization")
    st.plotly_chart(fig, use_container_width=True)

    # No
    st.subheader("Done regulatory clustering! Download your results using the button above.")



def run_run():
    thread = Thread(target=rcm_runner, args=())
    add_script_run_ctx(thread)
    thread.start()
    thread.join()


with c0:
    if methylation is not None and rna is not None and protein is not None:
        gene_id = st.selectbox('Gene identifier', list(set(methylation.columns) & set(rna.columns) & set(protein.columns)))
        st.button('Run RCM!', on_click=run_run)

    else:
        st.info(
            f"""
               Upload an DNA methylation dCpG file, a RNA file, and a protein file. 
            """
        )
        st.stop()