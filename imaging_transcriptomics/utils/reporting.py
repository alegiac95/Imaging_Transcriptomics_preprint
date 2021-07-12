import os
from pathlib import Path
from datetime import datetime
# Third party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from fpdf import FPDF

# TODO: update the info on the pdf


class PDF(FPDF):
    """Class to generate a PDF report for the imaging-transcriptomics script.
    """

    def header(self):
        """The header will contain always the title."""
        self.rect(10.0, 10.0, 190.0, 280.0)
        self.line(10.0, 50.0, 200.0, 50.0)
        self.set_font("Helvetica", "B", 14)
        self.cell(w=0, h=15.0, align="C",
                  txt="Virtual Histology Analysis Report", ln=True)

    def analysis_info(self, filename, date, filepath):
        """Info on the analysis performed. Information included are the name of
        the scan of the input, date of the analysis and the original path of the scan."""
        self.set_font("Courier", "", 10)
        self.cell(w=100.0, h=8.0, align="L", txt=f"  Scan Name: {filename}")
        self.cell(w=100.0, h=8.0, align="L", txt=f"  Date: {date}", ln=True)
        self.cell(w=100.0, h=10.0, align="L", txt=f"  File Path: {filepath}")

    def pls_regression(self, path_plots):
        """Include the plots of the pls components."""
        self.ln(20)
        self.set_font("Helvetica", "BU", 12)
        self.cell(w=0, h=10.0, align="L", txt="-PLS Regression")
        self.ln(10)
        self.image(Path(path_plots) / "variance.png", w=120.0)
        self.image(Path(path_plots) / "perc_variance.png", w=120.0)


def plot_variance(varexp, dim, sv_path):
    """Generate the plots to use in the pdf.
    The plots are sved as a png to eventually use in presentations and publications.

    INPUTS:
    - varexp: array of variance explained.
    - dim: optimal dimension selected to explain the variance we want.
    - sv_path: path where to save the plots.
    """
    # Plot variance explained --> styling of the plot needed
    save_path = Path(sv_path)
    perc_explained_var = 100 * np.cumsum(varexp)
    plt.plot(range(1, 16), 100 * np.cumsum(varexp),
             marker="o", color="sandybrown")
    plt.plot(dim, perc_explained_var[dim-1], 'o', color="red")
    plt.vlines(
        dim, perc_explained_var[0]-10, perc_explained_var[dim-1], colors="lightgrey", linestyles="dashed")
    plt.hlines(perc_explained_var[dim-1], 0, dim,
               colors="lightgrey", linestyles="dashed")
    plt.title("Cumulative variance explained by PLS components")
    plt.ylabel("Total explained variance (%)")
    plt.xlabel("Number of PLS components")
    plt.xlim(0, 15)
    plt.ylim(perc_explained_var[0]-10, 105)
    plt.grid(True)
    plt.savefig(save_path / "perc_variance.png", dpi=1200)
    plt.close()

    plt.bar(range(1, 16), 100 * varexp, color="sandybrown")
    for index, value in enumerate(varexp):
        plt.text(index + 0.5, 100 * value, "{:.1f}".format(100*value))
    plt.xlabel("PLS components")
    plt.ylabel("Variance (%)")
    plt.savefig(save_path / "variance.png", dpi=1200)
    plt.close()

    return


def reporting(scan_path,  varexp, dim, z, p_val, p_val_corr, pls, output_path=None):
    """Create the reporting for virtual histology.
    The repoting genereated is included in a folder named vh_{filename} and includes:
    - pdf with pls regression components and variance explained by the components.
    - individual csv files with z scores, ROI values, gene ids ans p values for each


    INPUTS:
    - scan_path: path of the scan used for the analysis. This is needed to generate the report folder name.
    - varexp:
    - output_path: path where the results are saved. If no path is provided the pathof the scan is used.
    """

    print("Creating report...")
    if output_path == None:
        output_path = Path(scan_path).parent
    else:
        output_path = Path(output_path)
    scan_path = Path(scan_path)

    scan_name = scan_path.name.split(".")[0]
    save_folder_name = f"vh_{scan_name}"
    matches = list(output_path.glob(f"{save_folder_name}*"))
    if matches:
        nr = len(matches)
        save_folder_name = f"{save_folder_name}_{nr}"
    out_directory = output_path / save_folder_name
    out_directory.mkdir()

    # Create the plots
    plot_variance(varexp, dim, out_directory)

    # Create the PDF report
    report_path = out_directory / "report.pdf"
    date = datetime.now().strftime("%d-%m-%Y")
    report = PDF(orientation="P", unit="mm", format="A4")
    report.add_page()
    report.analysis_info(filename=scan_path.name,
                         date=date, filepath=scan_path)
    report.pls_regression(path_plots=out_directory)
    report.output(report_path, 'F')

    # create the csv file(s) to save
    for i in z.keys():
        data = np.vstack((pls[i].reshape(1, 15633), z[i].reshape(
            1, 15633), p_val[i].reshape(1, 15633), p_val_corr[i].reshape(1, 15633))).T
        data = pd.DataFrame(data, columns=["Gene ID", "Z", "p", "p corrected"])
        data.to_csv(f"{out_directory}/PLS{i}.csv", index=False)


def print_table(var, p):
    print("+-----------+----------------+-------+")
    print("| Component | Cumulative var | p val |")
    print("|-----------|----------------|-------|")
    for i in range(var.size):
        print("|     {}     |      {:.3f}    | {} |".format(
            i+1, var[i], p[i]))
    print("+-----------+----------------+-------+")
    print("")
