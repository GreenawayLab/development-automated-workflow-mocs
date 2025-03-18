import os
import json
import pandas as pd
from pathlib import Path
import cagey 
from cagey import NmrPeak, NmrSpectrum
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import replace
from operator import attrgetter
from typing import Annotated
import typer

import nmrglue
import numpy as np
from rich.console import Console
from rich.table import Table


def get_spectrum(spectrum_dir: Path) -> NmrSpectrum:
    """Get NMR spectrum from the machine data directory.

    Parameters:
        spectrum_dir: Path to the directory containing the spectrum data.

    Returns:
        The NMR spectrum.
    """
    peaks = tuple(_pick_peaks(spectrum_dir))

    reference_peak_ppm = 1.99 # 1.99 previously
    reference_peak = _find_reference_peak(peaks, reference_peak_ppm)
    reference_shift = 1.96 - reference_peak.ppm
    acetonitrile_peaks = [1.96, 1.95, 1.94, 1.93]

    return NmrSpectrum(
        aldehyde_peaks=tuple(
            _get_aldehyde_peaks(peaks, reference_shift, acetonitrile_peaks)
        ),
        imine_peaks=tuple(
            _get_imine_peaks(peaks, reference_shift, acetonitrile_peaks)
        ),
    )

def _find_reference_peak(peaks: Sequence[NmrPeak], reference_peak_ppm: float) -> NmrPeak:
    """Find the reference peak closest to the specified ppm. Assigns the nearest peak 
    the value of the reference_peak_ppm if no exact match is found.
    
    Parameters:
        peaks: The sequence of detected NMR peaks.
        reference_peak_ppm: The target ppm value for the reference peak.

    Returns:
        The NmrPeak that is closest to the reference_peak_ppm.
    """
    possible_reference_peaks = list(
        filter(lambda peak: peak.has_ppm(reference_peak_ppm), peaks)
    )
    
    if not possible_reference_peaks:
        # Find the peak with the closest ppm to the reference
        closest_peak = min(peaks, key=lambda peak: abs(peak.ppm - reference_peak_ppm))
        # Assign it the reference ppm value
        closest_peak = replace(closest_peak, ppm=reference_peak_ppm)
        return closest_peak

    # If there's a match, return the peak with the highest amplitude
    return max(possible_reference_peaks, key=attrgetter("amplitude"))

def _pick_peaks(spectrum_dir: Path) -> Iterator[NmrPeak]:
    metadata, data = nmrglue.bruker.read_pdata(str(spectrum_dir))
    udic = nmrglue.bruker.guess_udic(metadata, data)
    unit_conversion = nmrglue.fileio.fileiobase.uc_from_udic(udic)
    for peak in nmrglue.peakpick.pick(data, pthres=1e4, nthres=None): #check threshold #1e4
        shift = unit_conversion.ppm(peak["X_AXIS"])
        amplitude = peak["VOL"]
        yield NmrPeak(shift, amplitude)

def _shift_peaks(
    peaks: Iterable[NmrPeak],
    shift: float,
) -> Iterator[NmrPeak]:
    for peak in peaks:
        yield replace(peak, ppm=peak.ppm + shift)


def _remove_peaks(
    peaks: Iterable[NmrPeak],
    to_remove: Sequence[float],
) -> Iterator[NmrPeak]:
    for peak in peaks:
        if not np.any(np.isclose(peak.ppm, to_remove, atol=0.03)):
            yield peak

def _get_aldehyde_peaks(
    peaks: Iterable[NmrPeak],
    reference_shift: float,
    solvent_peaks: Sequence[float],
) -> Iterator[NmrPeak]:
    peaks = filter(
        lambda peak: peak.in_range(9.1, 11.0), 
        peaks,
    )
    peaks = _shift_peaks(peaks, reference_shift)
    peaks = _remove_peaks(peaks, solvent_peaks)
    return filter(
        lambda peak: peak.amplitude > 0,
        peaks,
    )

def _get_imine_peaks(
    peaks: Iterable[NmrPeak],
    reference_shift: float,
    solvent_peaks: Sequence[float],
) -> Iterator[NmrPeak]:
    peaks = filter(
        lambda peak: peak.in_range(8.4, 8.8), #6.0, 9.0
        peaks,
    )
    peaks = _shift_peaks(peaks, reference_shift)
    peaks = _remove_peaks(peaks, solvent_peaks)
    return filter(
        lambda peak: peak.amplitude > 0,
        peaks,
    )

def calculate_aldehyde_percentage_per_peak(spectrum: NmrSpectrum) -> dict:
    """Calculate the percentage of each aldehyde peak compared to the highest imine peak.

    Parameters:
        spectrum: The NMR spectrum containing aldehyde and imine peaks.

    Returns:
        A dictionary with aldehyde peak details as keys and the percentage of their amplitude
        relative to the highest imine peak as values.
    """
    percentages = {}
    if not spectrum.imine_peaks:
        return percentages  # No imine peaks, return empty dict

    max_imine_peak = max(spectrum.imine_peaks, key=attrgetter("amplitude"))
    #print(max_imine_peak.amplitude)

    if max_imine_peak.amplitude == 0:
        return percentages  # No valid imine peak with amplitude

    for peak in spectrum.aldehyde_peaks:
        aldehyde_percentage = round((peak.amplitude / max_imine_peak.amplitude) * 100, 3)
        percentages[f"aldehyde_peak_{peak.ppm:.3f}"] = aldehyde_percentage

    return percentages



def add_aldehyde_comments(percentages: dict) -> dict:
    """Add comments to the aldehyde percentages based on the value.

    Parameters:
        percentages: A dictionary with aldehyde peak percentages.

    Returns:
        A dictionary with added comments based on the aldehyde percentage value.
    """
    result = {}
    for key, value in percentages.items():
        if value <= 6:
            comment = "aldehyde is below 5% and categorised as minimal aldehyde"
        elif 6 < value <= 20:
            comment = "aldehyde is between 5-20% and categorised as residual aldehyde"
        else:
            comment = "aldehyde is above 20% and categorised as significant aldehyde"
        
        result[key] = {
            "percentage": value,
            "aldehyde_comment": comment
        }

    return result

def process_nmr_file(title_file: Path) -> None:
    """Process a single NMR title file and generate a CSV file and JSON results."""
    title_file = Path(title_file)
    with open(title_file, 'r') as stream1:
        title_id = stream1.readline().strip()
    
    spectrum = get_spectrum(title_file.parent)

    aldehyde_data = {
        "id": [],
        "type": [],
        "ppm": [],
        "amplitude": [],
    }

    for id_, peak in enumerate(spectrum.aldehyde_peaks):
        aldehyde_data["id"].append(str(id_))
        aldehyde_data["type"].append("aldehyde")
        aldehyde_data["ppm"].append(str(peak.ppm))
        aldehyde_data["amplitude"].append(str(peak.amplitude))

    imine_data = {
        "id": [],
        "type": [],
        "ppm": [],
        "amplitude": [],
    }

    for id_, peak in enumerate(spectrum.imine_peaks):
        imine_data["id"].append(str(id_))
        imine_data["type"].append("imine")
        imine_data["ppm"].append(str(peak.ppm))
        imine_data["amplitude"].append(str(peak.amplitude))

    # Calculate the aldehyde percentage for each aldehyde peak compared to the highest imine peak
    aldehyde_percentages = calculate_aldehyde_percentage_per_peak(spectrum)
    
    # Add comments to the aldehyde percentages based on the calculated values
    aldehyde_results_with_comments = add_aldehyde_comments(aldehyde_percentages)
    
    # Save the results to a JSON file
    result_json_file = f'/home/abasford/projects/MOC_paper/metal_screen/nmr/{title_id}_result.json'
    with open(result_json_file, 'w') as json_file:
        json.dump(aldehyde_results_with_comments, json_file, indent=4)

    # Create DataFrames and save CSV file
    aldehyde_df = pd.DataFrame(aldehyde_data)
    imine_df = pd.DataFrame(imine_data)
    

    # Save the results to a JSON file
    result_json_file = f'/home/abasford/projects/MOC_paper/metal_screen/nmr/{title_id}_result.json'
    with open(result_json_file, 'w') as json_file:
        json.dump(aldehyde_results_with_comments, json_file, indent=4)

    result_nmr = pd.concat([aldehyde_df, imine_df], ignore_index=True)
    result_nmr.to_csv(f'/home/abasford/projects/MOC_paper/metal_screen/nmr/{title_id}.csv', index=False)

def process_directory(data_dir: Path) -> None:
    """Process all NMR files in the given directory structure.

    Parameters:
        data_dir: Path to the root directory containing multiple NMR data folders.
    """
    data_dir = Path(data_dir)

    # Walk through the directory structure
    for root, dirs, files in os.walk(data_dir):
        for file in files:
            if file == "title":
                title_file = Path(root) / file
                process_nmr_file(title_file)

def main(
    data_dir: Annotated[
        Path, typer.Argument(help="Path to the directory containing NMR data folders.")
    ]
) -> None:
    """Process all NMR data in the given directory."""
    process_directory(data_dir)

if __name__ == "__main__":
    typer.run(main)
