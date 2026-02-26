import uproot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import hashlib
import sys
import os
from datetime import datetime


class Waveforms:

    def __init__(self, filename="test5.root"):  # Initialize the class
        self.file = uproot.open(filename)  # Open the ROOT file
        self.events = self.file["tree"][0].arrays(
            library="pd"
        )  # Get the data from the TTree including events
        print(self.events.columns)
        self.TOPwv = self.file["tree"][8].arrays(
            library="pd"
        )  # Get the data from the TTree including waveforms
        print(self.TOPwv.columns)
        self.TOPrd = self.file["tree"][5].arrays(
            library="pd"
        )  # Get the data from the TTree including information on asic, carrier, and channel
        print(self.TOPrd.columns)
        self.propCycle = plt.rcParams["axes.prop_cycle"]
        self.colors = self.propCycle.by_key()[
            "color"
        ]  # Get the colors from the propCycle
        self.N = 50  # Number of events to plot
        # Get the first N unique event numbers
        unique_events = list(set(self.events.m_event))[: self.N]
        self.eventNumber = unique_events
        self.carriers = [
            0,
            1,
            2,
            3,
        ]  # List of carriers, there are only 4 on each boardstack
        self.asics = [0, 1, 2, 3]  # List of asics, there are only 4 on each carrier
        self.channels = 8  # Number of channels per asic

    def gettingWaveforms(self):  # Function to get the waveforms
        # Set up logging to file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_filename = f"waveform_analysis_log_{timestamp}.txt"

        # Create a custom print function that writes to both console and file
        log_file = open(log_filename, "w")
        original_stdout = sys.stdout

        class Logger:
            def __init__(self, console, file):
                self.console = console
                self.file = file

            def write(self, message):
                self.console.write(message)
                self.file.write(message)
                self.file.flush()  # Ensure immediate writing

            def flush(self):
                self.console.flush()
                self.file.flush()

        sys.stdout = Logger(original_stdout, log_file)

        print(f"Starting waveform analysis at {datetime.now()}")
        print(f"Logging output to: {log_filename}")
        print("=" * 50)

        plt.close("all")  # Rid of any leftovers
        fig, ax = plt.subplots(figsize=(30, 15))

        # Loop through all rows up to self.N to get data from the selected events
        count = 0
        wvHashes = set()
        all_channel_hits = []  # Storing the hits per chanel for CSV

        for i in range(len(self.events)):
            if count >= self.N:
                break

            ev = self.events.iloc[i].m_event

            # Skip if this event is not in our selected event numbers
            if ev not in self.eventNumber:
                continue

            # Get the corresponding data from all three TTrees using the same index i
            event_data = self.events.iloc[i]
            waveform_data = self.TOPwv.iloc[i]
            hardware_data = self.TOPrd.iloc[i]

            # Printing the row and event
            print(f"Row {i}, Event {ev}")

            # Extract hardware arrays from the hardware TTree
            def safe_extract_array(data, key):
                try:
                    val = data.get(key, [])
                    print(f"{key}: {val} (type: {type(val)})")
                    if hasattr(val, "tolist"):
                        return val.tolist()
                    elif hasattr(val, "__iter__"):
                        return list(val)
                    else:
                        return [val]
                except Exception as e:
                    print(f"Failed to extract {key}: {e}")
                    return []

            carriers = safe_extract_array(hardware_data, "TOPRawDigits.m_carrier")
            asics = safe_extract_array(hardware_data, "TOPRawDigits.m_asic")
            channels = safe_extract_array(hardware_data, "TOPRawDigits.m_channel")

            # Get waveform data arrays
            waveform_column = "TOPRawWaveforms.m_data"
            if waveform_column in self.TOPwv.columns:
                waveform_raw = waveform_data[waveform_column]
            else:
                waveform_raw = waveform_data.get(
                    "m_data", waveform_data.get("waveform", [])
                )

            print(f"Waveform raw type: {type(waveform_raw)}")
            print(
                f"Waveform raw shape: {getattr(waveform_raw, 'shape', 'No shape attr')}"
            )
            print(
                f"Waveform raw length: {len(waveform_raw) if hasattr(waveform_raw, '__len__') else 'No length'}"
            )

            if hasattr(waveform_raw, "__getitem__"):
                try:
                    print(
                        f"First element of raw waveform: shape/length = {getattr(waveform_raw[0], 'shape', len(waveform_raw[0]) if hasattr(waveform_raw[0], '__len__') else 'N/A')}"
                    )
                    print(
                        f"Second element of raw waveform: shape/length = {getattr(waveform_raw[1], 'shape', len(waveform_raw[1]) if hasattr(waveform_raw[1], '__len__') else 'N/A')}"
                    )
                except:
                    print("Could not access first/second element of raw waveform")
            print(
                f"Number of hardware entries: {len(carriers)} carriers, {len(asics)} asics, {len(channels)} channels"
            )

            # Convert waveforms - handle multiple channels properly
            try:
                if hasattr(waveform_raw, "tolist"):
                    # Convert awkward array to nested list structure
                    raw_list = waveform_raw.tolist()
                    print(
                        f"Raw list structure: {type(raw_list)}, length: {len(raw_list) if hasattr(raw_list, '__len__') else 'N/A'}"
                    )
                    if hasattr(raw_list, "__len__") and len(raw_list) > 0:
                        print(
                            f"First element: {type(raw_list[0])}, length: {len(raw_list[0]) if hasattr(raw_list[0], '__len__') else 'N/A'}"
                        )

                    # raw_list[0] and raw_list[1] are the two halves of all waveforms
                    # Need to concatenate them to get full 64-sample waveforms
                    if len(raw_list) >= 2:
                        # Concatenate the two halves
                        first_half = raw_list[0]
                        second_half = raw_list[1]
                        full_waveform = np.concatenate([first_half, second_half])
                        waveforms = [full_waveform]
                        print(
                            f"    Concatenated waveform: {len(first_half)} + {len(second_half)} = {len(full_waveform)} samples"
                        )
                    else:
                        # Only one part available
                        waveforms = [raw_list[0]] if len(raw_list) > 0 else []
                        print(
                            f"    Single waveform part: {len(raw_list[0]) if len(raw_list) > 0 else 0} samples"
                        )

                elif hasattr(waveform_raw, "__len__") and hasattr(
                    waveform_raw, "__getitem__"
                ):
                    # Handle direct access to waveform elements
                    waveforms = []
                    for j in range(len(waveform_raw)):
                        ch_waveform = waveform_raw[j]
                        if hasattr(ch_waveform, "tolist"):
                            waveforms.append(ch_waveform.tolist())
                        else:
                            waveforms.append(ch_waveform)
                        print(
                            f"    Channel {j} waveform: {len(waveforms[-1]) if hasattr(waveforms[-1], '__len__') else 'scalar'} samples"
                        )
                else:
                    waveforms = [waveform_raw]
                    print(
                        f"Direct waveform assignment: {len(waveform_raw) if hasattr(waveform_raw, '__len__') else 'scalar'} samples"
                    )

            except Exception as e:
                print(f"Waveform conversion failed: {e}")
                waveforms = []

            # First pass: collect all hits for this event and find the highest amplitude per channel
            num_channels = min(len(carriers), len(asics), len(channels), len(waveforms))
            print(f"Processing {num_channels} channels in this row")

            channel_hits = {}  # Disctionary for the current event

            for ch_idx in range(num_channels):
                carrier_val = int(carriers[ch_idx])
                asic_val = int(asics[ch_idx])
                channel_val = int(channels[ch_idx])
                waveform = waveforms[ch_idx]

                print(
                    f"Channel {ch_idx}: carrier={carrier_val}, asic={asic_val}, channel={channel_val}"
                )

                # Skip if this combination doesn't match what we want to plot
                if (
                    carrier_val not in self.carriers
                    or asic_val not in self.asics
                    or channel_val >= self.channels
                ):
                    continue

                # Waveform should now be a single concatenated list of 64 samples
                y = np.array(waveform) if waveform else np.array([])
                print(f"    Final waveform: {len(y)} samples")

                if (
                    len(y) == 0 or np.max(y) < 50
                ):  # Cut channels without peaks over 50 ADC counts
                    print(
                        f"    Skipping channel {carrier_val}-{asic_val}-{channel_val} due to low peak (max={np.max(y) if len(y) > 0 else 0:.1f})"
                    )
                    continue

                # Create simple x-axis starting from 0
                sampleNumber = len(y)
                x = np.arange(sampleNumber)

                channel_id = (carrier_val, asic_val, channel_val)
                max_amplitude = np.max(y)

                # Keep only the highest amplitude hit per channel for this event
                if (
                    channel_id not in channel_hits
                    or max_amplitude > channel_hits[channel_id]["max_amp"]
                ):
                    channel_hits[channel_id] = {
                        "x": x.copy(),
                        "y": y.copy(),
                        "carrier": carrier_val,
                        "asic": asic_val,
                        "channel": channel_val,
                        "max_amp": max_amplitude,
                        "event": ev,
                        "row": i,
                    }
                    print(
                        f"    New best hit for channel {carrier_val}-{asic_val}-{channel_val}: max={max_amplitude:.1f}"
                    )
                else:
                    print(
                        f"    Skipping lower amplitude hit for channel {carrier_val}-{asic_val}-{channel_val}: max={max_amplitude:.1f} < {channel_hits[channel_id]['max_amp']:.1f}"
                    )

            # Second pass: plot the best hit from each channel in this event
            for channel_id, hit_data in channel_hits.items():
                if count >= self.N:
                    break

                x = hit_data["x"]
                y = hit_data["y"]
                carrier_val = hit_data["carrier"]
                asic_val = hit_data["asic"]
                channel_val = hit_data["channel"]

                # Add minimal vertical separation between channels
                seq_ch = carrier_val * 32 + asic_val * 8 + channel_val
                threshold = 50  # Smaller separation
                y_offset = (
                    count * threshold
                )  # Use count instead of seq_ch for consistent spacing
                y = y + y_offset

                # Create a hash of the combined x,y data to delete repeated waveforms
                xHash = hashlib.md5(x.tobytes()).hexdigest()
                yHash = hashlib.md5(y.tobytes()).hexdigest()
                combinedHash = f"{xHash}_{yHash}"

                if combinedHash in wvHashes:
                    print("Skipping repeating waveform")
                    continue
                else:
                    wvHashes.add(combinedHash)
                    count += 1

                    # Store data for CSV export
                    all_channel_hits.append(
                        {
                            "waveform_number": count,
                            "event": hit_data["event"],
                            "row": hit_data["row"],
                            "carrier": carrier_val,
                            "asic": asic_val,
                            "channel": channel_val,
                            "seq_ch": seq_ch,
                            "max_amplitude": hit_data["max_amp"],
                            "waveform_length": len(x),
                        }
                    )

                    # Creating label
                    label = f"c{carrier_val} a{asic_val} ch{channel_val} ({seq_ch})"

                    try:
                        ax.plot(
                            x,
                            y,
                            linestyle="-",
                            linewidth=2,
                            color=self.colors[count % 10],
                            label=label,
                        )
                        print(
                            f"Plotted waveform {count}: event={hit_data['event']}, hw={carrier_val}-{asic_val}-{channel_val}, max_amp={hit_data['max_amp']:.1f}"
                        )

                    except Exception as e:
                        print(f"Error plotting waveform: {e}")

            # Break out of main loop if we've reached the limit
            if count >= self.N:
                break

        print(f"Total waveforms plotted: {count}")

        # Save channel hits as CSV
        if all_channel_hits:
            df_channel_hits = pd.DataFrame(all_channel_hits)
            df_channel_hits.to_csv("channel_hits.csv", index=False)
            print(f"Saved {len(all_channel_hits)} channel hits to channel_hits.csv")

        # Customizing the plot
        if count > 0:
            ax.legend(
                bbox_to_anchor=(1.02, 1),
                loc="upper left",
                fontsize=12,
                title="Channel Labels",
                title_fontsize=14,
                title_fontstyle="italic",
            )
        plt.title(
            f"{count} Channels with Hits Over 50 ADC Counts",
            fontsize=18,
            fontweight="bold",
        )
        plt.ylabel("Amplitude (with offset)", fontsize=14, fontweight="bold")
        plt.xlabel("Sample Number", fontsize=14, fontweight="bold")

        # Show grid for better visualization
        ax.grid(True, alpha=0.3)

        # Set reasonable axis limits
        if count > 0:
            plt.xlim(0, 64)  # Most waveforms are ~64 samples

        # Saving the plot to png
        plt.savefig("wvs-w-hits-over-50ADC.png", dpi=150)

        print("=" * 50)
        print(f"Analysis completed at {datetime.now()}")
        print(f"Log saved to: {log_filename}")

        # Restore original stdout and close log file
        sys.stdout = original_stdout
        log_file.close()

        print(f"All output has been saved to: {log_filename}")


def main():
    wv = Waveforms()
    pd.set_option("display.max_rows", None)  # So we can see all the rows
    pd.set_option(
        "display.max_columns", None
    )  # Uncomment this line if you want to see all columns

    # Print some basic info before plotting
    print(f"Total events in file: {len(wv.events)}")
    print(f"Unique events selected: {len(wv.eventNumber)}")
    print(f"Sample event numbers: {wv.eventNumber[:5]}")

    wv.gettingWaveforms()


if __name__ == "__main__":
    main()
