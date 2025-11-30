#!/usr/bin/env python3
"""
Fast mzML Peak Map Viewer using NiceGUI + Datashader + pyOpenMS

Designed to handle 50+ million peaks with smooth zooming and panning.
Uses datashader for server-side rendering of massive datasets.
"""

import io
import base64
import numpy as np
import pandas as pd
from pathlib import Path

# Datashader for fast rendering
import datashader as ds
import datashader.transfer_functions as tf
from colorcet import fire

# pyOpenMS for mzML loading
from pyopenms import MSExperiment, MzMLFile

# NiceGUI for the web interface
from nicegui import ui, app


class MzMLViewer:
    """High-performance mzML peak map viewer using datashader."""

    def __init__(self):
        self.exp = None
        self.df = None  # DataFrame with rt, mz, intensity
        self.current_file = None

        # View bounds (will be set after loading)
        self.rt_min = 0
        self.rt_max = 1
        self.mz_min = 0
        self.mz_max = 1

        # Current view (for zooming)
        self.view_rt_min = None
        self.view_rt_max = None
        self.view_mz_min = None
        self.view_mz_max = None

        # Image dimensions
        self.plot_width = 1200
        self.plot_height = 600

        # UI elements
        self.image_element = None
        self.status_label = None
        self.info_label = None
        self.rt_range_label = None
        self.mz_range_label = None

    def load_mzml(self, filepath: str) -> bool:
        """Load mzML file and extract peak data into a pandas DataFrame."""
        try:
            self.status_label.set_text(f"Loading {Path(filepath).name}...")
            ui.notify(f"Loading {filepath}...", type="info")

            # Load with pyOpenMS
            self.exp = MSExperiment()
            MzMLFile().load(filepath, self.exp)

            self.status_label.set_text("Extracting peaks...")

            # Pre-calculate total peaks for array allocation
            total_peaks = sum(spec.size() for spec in self.exp)

            if total_peaks == 0:
                ui.notify("No peaks found in file!", type="warning")
                return False

            # Pre-allocate arrays for speed
            rts = np.empty(total_peaks, dtype=np.float32)
            mzs = np.empty(total_peaks, dtype=np.float32)
            intensities = np.empty(total_peaks, dtype=np.float32)

            # Extract all peaks
            idx = 0
            for spec in self.exp:
                if spec.getMSLevel() != 1:
                    continue
                rt = spec.getRT()
                mz_array, int_array = spec.get_peaks()
                n = len(mz_array)
                if n > 0:
                    rts[idx:idx+n] = rt
                    mzs[idx:idx+n] = mz_array
                    intensities[idx:idx+n] = int_array
                    idx += n

            # Trim arrays to actual size (MS1 only)
            rts = rts[:idx]
            mzs = mzs[:idx]
            intensities = intensities[:idx]

            # Create DataFrame
            self.df = pd.DataFrame({
                'rt': rts,
                'mz': mzs,
                'intensity': intensities
            })

            # Log-transform intensity for better visualization
            self.df['log_intensity'] = np.log1p(self.df['intensity'])

            # Set bounds
            self.rt_min = float(self.df['rt'].min())
            self.rt_max = float(self.df['rt'].max())
            self.mz_min = float(self.df['mz'].min())
            self.mz_max = float(self.df['mz'].max())

            # Reset view to full range
            self.view_rt_min = self.rt_min
            self.view_rt_max = self.rt_max
            self.view_mz_min = self.mz_min
            self.view_mz_max = self.mz_max

            self.current_file = filepath

            # Update info
            self.info_label.set_text(
                f"Loaded: {Path(filepath).name} | "
                f"Spectra: {len(self.exp):,} | "
                f"Peaks: {len(self.df):,}"
            )
            self.status_label.set_text("Ready")
            ui.notify(f"Loaded {len(self.df):,} peaks", type="positive")

            return True

        except Exception as e:
            self.status_label.set_text(f"Error: {e}")
            ui.notify(f"Error loading file: {e}", type="negative")
            return False

    def render_image(self) -> str:
        """Render the current view using datashader and return base64 PNG."""
        if self.df is None or len(self.df) == 0:
            return ""

        # Filter data to current view
        mask = (
            (self.df['rt'] >= self.view_rt_min) &
            (self.df['rt'] <= self.view_rt_max) &
            (self.df['mz'] >= self.view_mz_min) &
            (self.df['mz'] <= self.view_mz_max)
        )
        view_df = self.df[mask]

        if len(view_df) == 0:
            return ""

        # Create datashader canvas
        canvas = ds.Canvas(
            plot_width=self.plot_width,
            plot_height=self.plot_height,
            x_range=(self.view_rt_min, self.view_rt_max),
            y_range=(self.view_mz_min, self.view_mz_max)
        )

        # Aggregate using mean of log intensity
        agg = canvas.points(view_df, 'rt', 'mz', ds.mean('log_intensity'))

        # Apply colormap (fire is great for MS data)
        img = tf.shade(agg, cmap=fire, how='linear')
        img = tf.set_background(img, 'black')

        # Convert to PNG bytes
        pil_img = img.to_pil()
        buffer = io.BytesIO()
        pil_img.save(buffer, format='PNG')
        buffer.seek(0)

        # Return base64 encoded
        return base64.b64encode(buffer.getvalue()).decode('utf-8')

    def update_plot(self):
        """Update the displayed plot."""
        if self.df is None:
            return

        self.status_label.set_text("Rendering...")

        img_data = self.render_image()
        if img_data:
            self.image_element.set_source(f"data:image/png;base64,{img_data}")

        # Update range labels
        self.rt_range_label.set_text(f"RT: {self.view_rt_min:.2f} - {self.view_rt_max:.2f} s")
        self.mz_range_label.set_text(f"m/z: {self.view_mz_min:.2f} - {self.view_mz_max:.2f}")

        self.status_label.set_text("Ready")

    def reset_view(self):
        """Reset to full view."""
        if self.df is None:
            return
        self.view_rt_min = self.rt_min
        self.view_rt_max = self.rt_max
        self.view_mz_min = self.mz_min
        self.view_mz_max = self.mz_max
        self.update_plot()

    def zoom_in(self, factor=0.5):
        """Zoom in by factor (0.5 = zoom to 50% of current range)."""
        if self.df is None:
            return
        rt_center = (self.view_rt_min + self.view_rt_max) / 2
        mz_center = (self.view_mz_min + self.view_mz_max) / 2
        rt_range = (self.view_rt_max - self.view_rt_min) * factor / 2
        mz_range = (self.view_mz_max - self.view_mz_min) * factor / 2

        self.view_rt_min = rt_center - rt_range
        self.view_rt_max = rt_center + rt_range
        self.view_mz_min = mz_center - mz_range
        self.view_mz_max = mz_center + mz_range
        self.update_plot()

    def zoom_out(self, factor=2.0):
        """Zoom out by factor."""
        if self.df is None:
            return
        rt_center = (self.view_rt_min + self.view_rt_max) / 2
        mz_center = (self.view_mz_min + self.view_mz_max) / 2
        rt_range = (self.view_rt_max - self.view_rt_min) * factor / 2
        mz_range = (self.view_mz_max - self.view_mz_min) * factor / 2

        # Clamp to data bounds
        self.view_rt_min = max(self.rt_min, rt_center - rt_range)
        self.view_rt_max = min(self.rt_max, rt_center + rt_range)
        self.view_mz_min = max(self.mz_min, mz_center - mz_range)
        self.view_mz_max = min(self.mz_max, mz_center + mz_range)
        self.update_plot()

    def pan(self, rt_frac=0, mz_frac=0):
        """Pan by fraction of current view range."""
        if self.df is None:
            return
        rt_shift = (self.view_rt_max - self.view_rt_min) * rt_frac
        mz_shift = (self.view_mz_max - self.view_mz_min) * mz_frac

        # Check bounds
        if self.view_rt_min + rt_shift < self.rt_min:
            rt_shift = self.rt_min - self.view_rt_min
        if self.view_rt_max + rt_shift > self.rt_max:
            rt_shift = self.rt_max - self.view_rt_max
        if self.view_mz_min + mz_shift < self.mz_min:
            mz_shift = self.mz_min - self.view_mz_min
        if self.view_mz_max + mz_shift > self.mz_max:
            mz_shift = self.mz_max - self.view_mz_max

        self.view_rt_min += rt_shift
        self.view_rt_max += rt_shift
        self.view_mz_min += mz_shift
        self.view_mz_max += mz_shift
        self.update_plot()

    def apply_custom_range(self, rt_min, rt_max, mz_min, mz_max):
        """Apply custom RT and m/z ranges."""
        if self.df is None:
            return
        try:
            self.view_rt_min = max(self.rt_min, float(rt_min))
            self.view_rt_max = min(self.rt_max, float(rt_max))
            self.view_mz_min = max(self.mz_min, float(mz_min))
            self.view_mz_max = min(self.mz_max, float(mz_max))
            self.update_plot()
        except ValueError:
            ui.notify("Invalid range values", type="warning")


def create_ui():
    """Create the NiceGUI interface."""
    viewer = MzMLViewer()

    # Dark theme for better contrast with fire colormap
    ui.dark_mode().enable()

    with ui.column().classes('w-full items-center p-4'):
        ui.label('mzML Peak Map Viewer').classes('text-3xl font-bold mb-4')
        ui.label('High-performance visualization with Datashader + pyOpenMS').classes('text-gray-400 mb-4')

        # File upload section
        with ui.row().classes('w-full justify-center mb-4'):
            async def handle_upload(e):
                content = e.content.read()
                # Save to temp file
                temp_path = Path('/tmp') / e.name
                temp_path.write_bytes(content)
                if viewer.load_mzml(str(temp_path)):
                    viewer.update_plot()

            ui.upload(
                label='Upload mzML file',
                on_upload=handle_upload,
                auto_upload=True
            ).props('accept=.mzML,.mzml').classes('max-w-lg')

            # Or load from path
            file_input = ui.input(
                label='Or enter file path',
                placeholder='/path/to/file.mzML'
            ).classes('w-96')

            async def load_from_path():
                path = file_input.value
                if path and Path(path).exists():
                    if viewer.load_mzml(path):
                        viewer.update_plot()
                else:
                    ui.notify("File not found", type="warning")

            ui.button('Load', on_click=load_from_path).props('color=primary')

        # Info bar
        with ui.row().classes('w-full justify-center gap-8 mb-2'):
            viewer.info_label = ui.label('No file loaded').classes('text-gray-400')
            viewer.status_label = ui.label('Ready').classes('text-green-400')

        # Range display
        with ui.row().classes('w-full justify-center gap-8 mb-2'):
            viewer.rt_range_label = ui.label('RT: -- - -- s').classes('text-blue-300')
            viewer.mz_range_label = ui.label('m/z: -- - --').classes('text-blue-300')

        # The main plot image
        with ui.card().classes('p-0'):
            viewer.image_element = ui.image().classes('w-full').style(
                f'width: {viewer.plot_width}px; height: {viewer.plot_height}px; background: black;'
            )

        # Navigation controls
        with ui.row().classes('justify-center gap-2 mt-4'):
            ui.button('Reset View', on_click=viewer.reset_view).props('color=secondary')
            ui.button('Zoom In', on_click=lambda: viewer.zoom_in(0.5)).props('color=primary')
            ui.button('Zoom Out', on_click=lambda: viewer.zoom_out(2.0)).props('color=primary')

        with ui.row().classes('justify-center gap-2 mt-2'):
            ui.button('← Pan Left', on_click=lambda: viewer.pan(rt_frac=-0.25)).props('color=accent')
            ui.button('→ Pan Right', on_click=lambda: viewer.pan(rt_frac=0.25)).props('color=accent')
            ui.button('↑ Pan Up', on_click=lambda: viewer.pan(mz_frac=0.25)).props('color=accent')
            ui.button('↓ Pan Down', on_click=lambda: viewer.pan(mz_frac=-0.25)).props('color=accent')

        # Custom range inputs
        with ui.expansion('Custom Range', icon='tune').classes('w-full max-w-4xl mt-4'):
            with ui.row().classes('w-full gap-4 items-end'):
                rt_min_input = ui.number(label='RT Min (s)', value=0, format='%.2f')
                rt_max_input = ui.number(label='RT Max (s)', value=1000, format='%.2f')
                mz_min_input = ui.number(label='m/z Min', value=0, format='%.2f')
                mz_max_input = ui.number(label='m/z Max', value=2000, format='%.2f')

                def apply_range():
                    viewer.apply_custom_range(
                        rt_min_input.value, rt_max_input.value,
                        mz_min_input.value, mz_max_input.value
                    )

                ui.button('Apply Range', on_click=apply_range).props('color=primary')

        # Keyboard shortcuts info
        with ui.expansion('Keyboard Shortcuts', icon='keyboard').classes('w-full max-w-4xl mt-2'):
            ui.markdown('''
            | Key | Action |
            |-----|--------|
            | `+` / `=` | Zoom In |
            | `-` | Zoom Out |
            | `Arrow Keys` | Pan |
            | `Home` | Reset View |
            ''')

        # Add keyboard handlers
        ui.keyboard(
            on_key=lambda e: (
                viewer.zoom_in() if e.key in ['+', '='] and e.action.keydown else
                viewer.zoom_out() if e.key == '-' and e.action.keydown else
                viewer.pan(rt_frac=-0.1) if e.key.arrow_left and e.action.keydown else
                viewer.pan(rt_frac=0.1) if e.key.arrow_right and e.action.keydown else
                viewer.pan(mz_frac=0.1) if e.key.arrow_up and e.action.keydown else
                viewer.pan(mz_frac=-0.1) if e.key.arrow_down and e.action.keydown else
                viewer.reset_view() if e.key == 'Home' and e.action.keydown else
                None
            )
        )


# Create the UI
create_ui()

# Run the app
if __name__ in {"__main__", "__mp_main__"}:
    ui.run(
        title='mzML Peak Map Viewer',
        host='0.0.0.0',
        port=8080,
        reload=False,
        show=False
    )
