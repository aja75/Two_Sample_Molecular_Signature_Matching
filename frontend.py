import PySimpleGUI as sg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from datetime import datetime

import MSM

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib

matplotlib.use('TkAgg')

version = '3.28.25'

font = 'Helvetica'

sg.theme('DefaultNoMoreNagging')

current_canvas = None  # Tracks the current canvas image

# Define the layout of the title frame
title_font_size = 16
subtitle_font_size = 12

title_frame = [
    [sg.Text('Two Sample Molecular Signature Matching', font=(font, title_font_size))],
    [sg.Text('Developed by: Anthony Aceto', font=(font, subtitle_font_size))],
    [sg.Text(f'Version: {version}', font=(font, subtitle_font_size))]
]

instruction_font_size = 10

data_input_frame = [
    [sg.Text('Control Sample: ', font=(font, subtitle_font_size)),
     sg.Input(font=(font, subtitle_font_size), key='-CTRL_INPUT-'), sg.FileBrowse(font=(font, subtitle_font_size))],

    [sg.Text('Test Sample: ', font=(font, subtitle_font_size)),
     sg.Input(font=(font, subtitle_font_size), key='-TEST_INPUT-'), sg.FileBrowse(font=(font, subtitle_font_size))],

    [sg.Text('Control DEG Format: ', font=(font, subtitle_font_size)),
     sg.Radio('log2(fold change)', group_id=1, default=True, key='-CTRL_FC-'),
     sg.Radio('Z-Score', group_id=1, key='-CTRL_Z-'),
     sg.Text('Threshold: ', font=(font, subtitle_font_size)),
     sg.Input('1', size=(4, 1), key='-CTRL_FMT_THRESHOLD-')],

    [sg.Text('Test DEG Format: ', font=(font, subtitle_font_size)),
     sg.Radio('log2(fold change)', group_id=2, default=True, key='-TEST_FC-'),
     sg.Radio('Z-Score', group_id=2, key='-TEST_Z-'),
     sg.Text('Threshold: ', font=(font, subtitle_font_size)),
     sg.Input('1', size=(4, 1), key='-TEST_FMT_THRESHOLD-')],

    [sg.Text('Gene Column Name: ', font=(font, subtitle_font_size)), sg.Input('Gene', key='-GENE_COL_NAME-')],

    [sg.Text('Number of Permutations: ', font=(font, subtitle_font_size)), sg.Input('100', key='-NUM_PERT-')],

    [sg.Text('p-value Significance Threshold: ', font=(font, subtitle_font_size)),
     sg.Radio('a<=0.01', group_id=3, key='-0.01-'),
     sg.Radio('a<=0.05', group_id=3, default=True, key='-0.05-'),
     sg.Radio('a<=0.10', group_id=3, key='-0.10-'),
     sg.Radio('a<=0.25', group_id=3, key='-0.25-')],

]

result_frame = [
    [sg.Canvas(key='-CANVAS-', size=(400, 400)),  # Canvas on the left
     sg.Column([  # Permutation Statistics Table
         [sg.Text('Permutation Statistics', font=(font, subtitle_font_size), visible=True, key='-TABLE_TXT-')],
         [sg.Table(
             values=[],
             headings=['Metric', 'Value'],
             auto_size_columns=False,
             justification='left',
             num_rows=2,
             key='-STATS_TABLE-',
             enable_events=True,
             row_height=25,
             visible=True,
             hide_vertical_scroll=True,
             col_widths=[10, 15]
         )],
         [sg.Text('')],
         [sg.Text('DEG Binary Mapping', font=(font, subtitle_font_size), visible=True, key='-BINARY_TXT-')],
         [sg.Table(
                     values=[],
                     headings=['MAP', 'Control', 'Test', 'Mean', 'Std Dev'],
                     auto_size_columns=False,
                     justification='left',
                     num_rows=3,
                     key='-BINARY_TABLE-',
                     enable_events=True,
                     row_height=25,
                     visible=True,
                     hide_vertical_scroll=True,
                     col_widths=[5, 10, 10, 10, 10]
                 )],
         [sg.Button('Export', key='-EXPORT-')]
     ],
         visible=False, key='-TABLE_COL-')
     ],

]

main_layout = [
    [sg.Frame('', title_frame, element_justification='left', border_width=0)],
    [sg.HorizontalSeparator()],
    [sg.Frame('', data_input_frame, element_justification='right', border_width=0), sg.VerticalSeparator(),
     sg.Frame('', result_frame, element_justification='right', border_width=0)],
    [sg.Button("Submit"), sg.Button("Exit")]
]

window = sg.Window("Two Sample Molecular Signature Matching",
                   main_layout,
                   size=(1500, 600),
                   resizable=True,
                   finalize=True)


def draw_figure(canvas, figure):
    global current_canvas
    if current_canvas is not None:
        current_canvas.get_tk_widget().pack_forget()
        current_canvas = None

    current_canvas = FigureCanvasTkAgg(figure, canvas)
    current_canvas.draw()
    current_canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
    return current_canvas


def plot_outcome(plot_metrics, table_data, binary_data):
    fig = matplotlib.figure.Figure(figsize=(4, 4), dpi=100)
    ax = fig.add_subplot(111)

    ax.axhline(y=0, color='grey', linestyle='--', linewidth=1, alpha=0.5)

    ax.set_ylim(min(plot_metrics['Values'] - 2), max(plot_metrics['Values']) + 2)

    bars = ax.bar(plot_metrics['Metrics'],
                  plot_metrics['Values'],
                  color=['blue', 'green'])

    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2,  # x-coordinate (center of bar)
                height,  # y-coordinate (top of the bar)
                f'{height:.2f}%',  # Text (formatted value)
                ha='center',  # Horizontal alignment
                va='bottom',  # Vertical alignment
                fontsize=10)  # Font size

    ax.set_title("Metrics Overview")
    ax.set_ylabel("Percentage (%)")
    ax.set_xlabel('')

    global current_canvas
    current_canvas = draw_figure(window['-CANVAS-'].TKCanvas, fig)

    window['-TABLE_COL-'].update(visible=True)
    window['-STATS_TABLE-'].update(values=table_data)
    window['-BINARY_TABLE-'].update(values=binary_data)


# Main Event Loop
while True:
    event, values = window.read()

    # Close the window
    if event == sg.WIN_CLOSED or event == "Exit":
        break

    elif event == 'Submit':
        DEG_profiles_path = [str(values['-CTRL_INPUT-']), str(values['-TEST_INPUT-'])]
        gene_index_col = str(values['-GENE_COL_NAME-'])
        n_permutations = int(values['-NUM_PERT-'])

        sig_thresholds = ['-0.01-', '-0.05-', '-0.10-', '-0.25-']
        significance_threshold = float([key for key in sig_thresholds if values.get(key)][0].replace("-", ""))

        control_formats = ['-CTRL_FC-', '-CTRL_Z-']
        control_format = str([key for key in control_formats if values.get(key)][0].replace("-", ""))
        control_format = 'log2FC' if "FC" in control_format else 'Z-Score'

        control_threshold = float(values['-CTRL_FMT_THRESHOLD-'])

        test_formats = ['-TEST_FC-', '-TEST_Z-']
        test_format = str([key for key in test_formats if values.get(key)][0].replace("-", ""))
        test_format = 'log2FC' if "FC" in test_format else 'Z-Score'

        test_threshold = float(values['-TEST_FMT_THRESHOLD-'])

        DEG_formats = [[control_format, control_threshold], [test_format, test_threshold]]

        try:
            DEG_Overlap, observed_concordance, p_value, map_data = MSM.main(DEG_profiles_path=DEG_profiles_path,
                                                                            DEG_formats=DEG_formats,
                                                                            gene_index_col=gene_index_col,
                                                                            n_permutations=n_permutations,
                                                                            test_column_index=1)

            plot_metrics = pd.DataFrame({
                'Metrics': ['DEG Overlap', 'Concordance'],
                'Values': [DEG_Overlap * 100, observed_concordance * 100],
            })

            table_data = [['p-value', f'{p_value:.4e}'],
                          ['Significance',
                           f' {"Not Significant" if p_value > significance_threshold else "Significant"}']]

            map_data = map_data.reindex([0, 1, -1])

            binary_data = [['0', f'{map_data.iloc[0, 0]:.2f}%', f'{map_data.iloc[0, 1]:.2f}%',
                            f'{map_data.iloc[0, 2]:.2f}%', f'{map_data.iloc[0, 3]:.2f}%'],

                           ['1', f'{map_data.iloc[1, 0]:.2f}%', f'{map_data.iloc[1, 1]:.2f}%',
                            f'{map_data.iloc[1, 2]:.2f}%', f'{map_data.iloc[1, 3]:.2f}%'],

                           ['-1', f'{map_data.iloc[2, 0]:.2f}%', f'{map_data.iloc[2, 1]:.2f}%',
                            f'{map_data.iloc[2, 2]:.2f}%', f'{map_data.iloc[2, 3]:.2f}%']
                           ]

            plot_outcome(plot_metrics, table_data, binary_data)

        except Exception as e:
            sg.popup_ok(e)

        else:
            pass

    elif event == '-EXPORT-':
        output_directory = sg.popup_get_folder('Select Output Directory')

        output_data = pd.DataFrame({
            'Index': ['Control Source', 'Control DEG Format', 'Control DEG Threshold',
                      'Test DEG Format', 'Test DEG Threshold', 'Test Source',
                      'DEG Overlap', 'Concordance', 'n permutations',
                      'p-value', 'Significance Threshold', 'Significance'],
            'Value': [DEG_profiles_path[0], DEG_formats[0][0], DEG_formats[0][1],
                      DEG_formats[1][0], DEG_formats[1][1], DEG_profiles_path[1],
                      DEG_Overlap, observed_concordance, n_permutations, p_value,
                      f'a<={significance_threshold}',
                      f'{"Not Significant" if p_value > significance_threshold else "Significant"}']
        })

        output_data.to_csv(f'{output_directory}/{datetime.now().strftime("%M_%S")}_MSM_data.csv', index=False)


# Close the window when done
window.close()
