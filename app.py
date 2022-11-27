import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sircle import SciRCM
from matplotlib import rcParams
from threading import Thread
from sciviso import Sankeyplot

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

#st.set_page_config(page_icon="✂️", page_title="SiRCle")

# st.image(
#     "https://emojipedia-us.s3.dualstack.us-west-1.amazonaws.com/thumbs/240/apple/285/scissors_2702-fe0f.png",
#     width=100,
# )

st.title("Signature Regulatory Clustering Model (Beta mode)")

background_method = st.text_input('Background Method', 'P&R')
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

    download_button(
        df,
        "SiRCle.csv",
        "Download SiRCle clustering to CSV",
    )

    plot_df = df[df['Regulation_Grouping_2'] != 'None']
    sk = Sankeyplot(plot_df)
    fig = sk.plot(columns=['Methylation', 'RNA', 'Protein', 'Regulation_Grouping_2',
                           'Regulation_Grouping_2_filtered'], colour_col='Protein')
    st.plotly_chart(fig)
    plot_bar = False
    if plot_bar:
        # # figure size in inches
        rcParams['figure.figsize'] = 3, 2
        sns.set(rc={'figure.figsize': (3, 2)}, style='ticks')

        colour_map = {'MDS': '#d8419b', 'MDS_TMDE': '#e585c0', 'MDS_ncRNA': '#d880b4',
                      'MDE': '#6aaf44', 'MDE_TMDS': '#0e8e6d', 'MDE_ncRNA': '#9edb77',
                      'TMDE': '#fe2323', 'TMDS': '#2952ff',
                      'TPDE': '#e68e25', 'TPDE_TMDS': '#844c0f',
                      'TPDS': '#462d76', 'TPDS_TMDE': '#9b29b7'}

        sns.set_palette("Greys_r")
        rcm_labels = ["MDS", "MDS_TMDE", "MDE", "MDE_TMDS", "TMDE", "TMDS", "TPDE", "TPDE_TMDS", "TPDS", "TPDS_TMDE"]
        colours = [colour_map[c] for c in rcm_labels]

        g = sns.catplot(data=df, x='Regulation_Grouping_2', kind="count", order=rcm_labels,
                        palette=sns.color_palette(colours),
                        height=4)
        plt.xticks(rotation=45, ha='right')
        plt.title(f'Regulation_Grouping_2')
        st.pyplot(plt.gcf())

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