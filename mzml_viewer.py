#!/usr/bin/env python3
"""
Fast mzML Peak Map Viewer using NiceGUI + Datashader + pyOpenMS

Designed to handle 50+ million peaks with smooth zooming and panning.
Uses datashader for server-side rendering of massive datasets.
Supports FeatureMap overlay with centroids, bounding boxes, and convex hulls.
"""

import io
import base64
import math
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Optional

# Datashader for fast rendering
import datashader as ds
import datashader.transfer_functions as tf
from colorcet import fire

# PIL for drawing overlays and axes
from PIL import Image, ImageDraw, ImageFont

# pyOpenMS for mzML and featureXML loading
from pyopenms import MSExperiment, MzMLFile, FeatureMap, FeatureXMLFile

# NiceGUI for the web interface
from nicegui import ui, app


def calculate_nice_ticks(vmin: float, vmax: float, num_ticks: int = 6) -> List[float]:
    """Calculate nice round tick values for an axis."""
    if vmin >= vmax:
        return [vmin]

    range_val = vmax - vmin
    rough_step = range_val / (num_ticks - 1)

    # Find the order of magnitude
    mag = math.floor(math.log10(rough_step))
    pow10 = 10 ** mag

    # Normalize step to 1-10 range
    norm_step = rough_step / pow10

    # Choose a nice step value
    if norm_step < 1.5:
        nice_step = 1
    elif norm_step < 3:
        nice_step = 2
    elif norm_step < 7:
        nice_step = 5
    else:
        nice_step = 10

    step = nice_step * pow10

    # Calculate tick positions
    first_tick = math.ceil(vmin / step) * step
    ticks = []
    tick = first_tick
    while tick <= vmax + step * 0.001:  # Small tolerance for floating point
        ticks.append(tick)
        tick += step

    return ticks


def format_tick_label(value: float, range_val: float) -> str:
    """Format a tick label based on the value and range."""
    if range_val >= 1000:
        if abs(value) >= 1000:
            return f"{value:.0f}"
        return f"{value:.1f}"
    elif range_val >= 10:
        return f"{value:.1f}"
    elif range_val >= 1:
        return f"{value:.2f}"
    else:
        return f"{value:.3f}"


class MzMLViewer:
    """High-performance mzML peak map viewer using datashader with feature overlay."""

    def __init__(self):
        self.exp = None
        self.df = None  # DataFrame with rt, mz, intensity
        self.current_file = None

        # FeatureMap data
        self.feature_map = None
        self.features_file = None

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

        # Image dimensions (plot area only)
        self.plot_width = 1100
        self.plot_height = 550

        # Axis margins
        self.margin_left = 80     # Space for y-axis labels
        self.margin_right = 20
        self.margin_top = 20
        self.margin_bottom = 50   # Space for x-axis labels

        # Total canvas size
        self.canvas_width = self.plot_width + self.margin_left + self.margin_right
        self.canvas_height = self.plot_height + self.margin_top + self.margin_bottom

        # Feature display options
        self.show_centroids = True
        self.show_bounding_boxes = True
        self.show_convex_hulls = True

        # Feature colors (RGBA)
        self.centroid_color = (0, 255, 100, 255)  # Green
        self.bbox_color = (255, 255, 0, 200)  # Yellow
        self.hull_color = (0, 200, 255, 150)  # Cyan

        # Axis colors
        self.axis_color = (200, 200, 200, 255)  # Light gray
        self.tick_color = (180, 180, 180, 255)
        self.label_color = (220, 220, 220, 255)
        self.grid_color = (60, 60, 60, 255)  # Subtle grid

        # UI elements
        self.image_element = None
        self.status_label = None
        self.info_label = None
        self.feature_info_label = None
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

    def load_featuremap(self, filepath: str) -> bool:
        """Load featureXML file."""
        try:
            self.status_label.set_text(f"Loading features from {Path(filepath).name}...")
            ui.notify(f"Loading {filepath}...", type="info")

            self.feature_map = FeatureMap()
            FeatureXMLFile().load(filepath, self.feature_map)

            self.features_file = filepath

            n_features = self.feature_map.size()
            self.feature_info_label.set_text(f"Features: {n_features:,}")
            self.status_label.set_text("Ready")
            ui.notify(f"Loaded {n_features:,} features", type="positive")

            return True

        except Exception as e:
            self.status_label.set_text(f"Error: {e}")
            ui.notify(f"Error loading features: {e}", type="negative")
            return False

    def clear_features(self):
        """Clear loaded feature map."""
        self.feature_map = None
        self.features_file = None
        self.feature_info_label.set_text("Features: None")
        ui.notify("Features cleared", type="info")

    def _data_to_plot_pixel(self, rt: float, mz: float) -> Tuple[int, int]:
        """Convert RT/m/z coordinates to pixel coordinates within the plot area."""
        rt_range = self.view_rt_max - self.view_rt_min
        mz_range = self.view_mz_max - self.view_mz_min

        if rt_range == 0 or mz_range == 0:
            return (0, 0)

        x = int((rt - self.view_rt_min) / rt_range * self.plot_width)
        # Invert Y axis (higher m/z at top)
        y = int((1 - (mz - self.view_mz_min) / mz_range) * self.plot_height)

        return (x, y)

    def _plot_to_canvas_pixel(self, plot_x: int, plot_y: int) -> Tuple[int, int]:
        """Convert plot pixel coordinates to canvas pixel coordinates."""
        return (plot_x + self.margin_left, plot_y + self.margin_top)

    def _coord_to_pixel(self, rt: float, mz: float) -> Tuple[int, int]:
        """Convert RT/m/z coordinates to canvas pixel coordinates."""
        plot_x, plot_y = self._data_to_plot_pixel(rt, mz)
        return self._plot_to_canvas_pixel(plot_x, plot_y)

    def _is_in_view(self, rt: float, mz: float) -> bool:
        """Check if a point is within the current view."""
        return (self.view_rt_min <= rt <= self.view_rt_max and
                self.view_mz_min <= mz <= self.view_mz_max)

    def _feature_intersects_view(self, rt_min: float, rt_max: float,
                                  mz_min: float, mz_max: float) -> bool:
        """Check if a feature's bounding box intersects the current view."""
        return not (rt_max < self.view_rt_min or rt_min > self.view_rt_max or
                    mz_max < self.view_mz_min or mz_min > self.view_mz_max)

    def _draw_features_on_plot(self, img: Image.Image) -> Image.Image:
        """Draw feature overlays on the plot image (before adding axes)."""
        if self.feature_map is None or self.feature_map.size() == 0:
            return img

        # Convert to RGBA for transparency
        img = img.convert('RGBA')
        overlay = Image.new('RGBA', img.size, (0, 0, 0, 0))
        draw = ImageDraw.Draw(overlay)

        features_drawn = 0
        max_features = 10000  # Limit for performance

        for feature in self.feature_map:
            if features_drawn >= max_features:
                break

            # Get feature properties
            rt = feature.getRT()
            mz = feature.getMZ()

            # Get bounding box from convex hulls
            hulls = feature.getConvexHulls()
            if hulls:
                all_points = []
                for hull in hulls:
                    points = hull.getHullPoints()
                    all_points.extend([(p[0], p[1]) for p in points])

                if all_points:
                    rt_coords = [p[0] for p in all_points]
                    mz_coords = [p[1] for p in all_points]
                    feat_rt_min, feat_rt_max = min(rt_coords), max(rt_coords)
                    feat_mz_min, feat_mz_max = min(mz_coords), max(mz_coords)
                else:
                    feat_rt_min, feat_rt_max = rt - 1, rt + 1
                    feat_mz_min, feat_mz_max = mz - 0.5, mz + 0.5
            else:
                feat_rt_min, feat_rt_max = rt - 1, rt + 1
                feat_mz_min, feat_mz_max = mz - 0.5, mz + 0.5

            # Skip if not in view
            if not self._feature_intersects_view(feat_rt_min, feat_rt_max,
                                                  feat_mz_min, feat_mz_max):
                continue

            features_drawn += 1

            # Draw convex hulls (using plot coordinates, not canvas)
            if self.show_convex_hulls and hulls:
                for hull in hulls:
                    points = hull.getHullPoints()
                    if len(points) >= 3:
                        pixel_points = [self._data_to_plot_pixel(p[0], p[1]) for p in points]
                        pixel_points.append(pixel_points[0])
                        draw.polygon(pixel_points, outline=self.hull_color,
                                    fill=(*self.hull_color[:3], 50))

            # Draw bounding box
            if self.show_bounding_boxes:
                top_left = self._data_to_plot_pixel(feat_rt_min, feat_mz_max)
                bottom_right = self._data_to_plot_pixel(feat_rt_max, feat_mz_min)
                draw.rectangle([top_left, bottom_right], outline=self.bbox_color, width=1)

            # Draw centroid
            if self.show_centroids:
                cx, cy = self._data_to_plot_pixel(rt, mz)
                r = 3
                draw.ellipse([cx-r, cy-r, cx+r, cy+r], fill=self.centroid_color,
                            outline=(255, 255, 255, 255))

        img = Image.alpha_composite(img, overlay)
        return img

    def _draw_axes(self, canvas: Image.Image) -> Image.Image:
        """Draw axes, tick marks, labels, and axis titles on the canvas."""
        draw = ImageDraw.Draw(canvas)

        # Try to load a font, fall back to default
        try:
            font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 12)
            title_font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 14)
        except:
            try:
                font = ImageFont.truetype("/usr/share/fonts/TTF/DejaVuSans.ttf", 12)
                title_font = ImageFont.truetype("/usr/share/fonts/TTF/DejaVuSans.ttf", 14)
            except:
                font = ImageFont.load_default()
                title_font = font

        # Plot area boundaries on canvas
        plot_left = self.margin_left
        plot_right = self.margin_left + self.plot_width
        plot_top = self.margin_top
        plot_bottom = self.margin_top + self.plot_height

        # Draw plot border
        draw.rectangle(
            [plot_left, plot_top, plot_right, plot_bottom],
            outline=self.axis_color,
            width=1
        )

        # --- X-axis (RT) ---
        rt_ticks = calculate_nice_ticks(self.view_rt_min, self.view_rt_max, num_ticks=8)
        rt_range = self.view_rt_max - self.view_rt_min

        for tick_val in rt_ticks:
            if self.view_rt_min <= tick_val <= self.view_rt_max:
                # Calculate x position
                x_frac = (tick_val - self.view_rt_min) / rt_range
                x = plot_left + int(x_frac * self.plot_width)

                # Draw tick mark
                draw.line([(x, plot_bottom), (x, plot_bottom + 5)], fill=self.tick_color, width=1)

                # Draw grid line (subtle)
                draw.line([(x, plot_top), (x, plot_bottom)], fill=self.grid_color, width=1)

                # Draw label
                label = format_tick_label(tick_val, rt_range)
                bbox = draw.textbbox((0, 0), label, font=font)
                label_width = bbox[2] - bbox[0]
                draw.text((x - label_width // 2, plot_bottom + 8), label, fill=self.label_color, font=font)

        # X-axis title
        x_title = "RT (s)"
        bbox = draw.textbbox((0, 0), x_title, font=title_font)
        title_width = bbox[2] - bbox[0]
        draw.text(
            (plot_left + self.plot_width // 2 - title_width // 2, plot_bottom + 28),
            x_title, fill=self.label_color, font=title_font
        )

        # --- Y-axis (m/z) ---
        mz_ticks = calculate_nice_ticks(self.view_mz_min, self.view_mz_max, num_ticks=8)
        mz_range = self.view_mz_max - self.view_mz_min

        for tick_val in mz_ticks:
            if self.view_mz_min <= tick_val <= self.view_mz_max:
                # Calculate y position (inverted - higher m/z at top)
                y_frac = 1 - (tick_val - self.view_mz_min) / mz_range
                y = plot_top + int(y_frac * self.plot_height)

                # Draw tick mark
                draw.line([(plot_left - 5, y), (plot_left, y)], fill=self.tick_color, width=1)

                # Draw grid line (subtle)
                draw.line([(plot_left, y), (plot_right, y)], fill=self.grid_color, width=1)

                # Draw label
                label = format_tick_label(tick_val, mz_range)
                bbox = draw.textbbox((0, 0), label, font=font)
                label_width = bbox[2] - bbox[0]
                label_height = bbox[3] - bbox[1]
                draw.text((plot_left - label_width - 10, y - label_height // 2), label, fill=self.label_color, font=font)

        # Y-axis title (rotated text - we'll draw it vertically)
        y_title = "m/z"
        bbox = draw.textbbox((0, 0), y_title, font=title_font)
        title_height = bbox[3] - bbox[1]

        # Create a small image for rotated text
        txt_img = Image.new('RGBA', (100, 30), (0, 0, 0, 0))
        txt_draw = ImageDraw.Draw(txt_img)
        txt_draw.text((0, 0), y_title, fill=self.label_color, font=title_font)
        txt_img = txt_img.rotate(90, expand=True)

        # Paste rotated text
        y_title_x = 5
        y_title_y = plot_top + self.plot_height // 2 - txt_img.height // 2
        canvas.paste(txt_img, (y_title_x, y_title_y), txt_img)

        return canvas

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

        # Create datashader canvas for the plot area only
        ds_canvas = ds.Canvas(
            plot_width=self.plot_width,
            plot_height=self.plot_height,
            x_range=(self.view_rt_min, self.view_rt_max),
            y_range=(self.view_mz_min, self.view_mz_max)
        )

        # Aggregate using mean of log intensity
        agg = ds_canvas.points(view_df, 'rt', 'mz', ds.mean('log_intensity'))

        # Apply colormap
        img = tf.shade(agg, cmap=fire, how='linear')
        img = tf.set_background(img, 'black')

        # Convert to PIL Image
        plot_img = img.to_pil()

        # Draw feature overlays on the plot image
        if self.feature_map is not None:
            plot_img = self._draw_features_on_plot(plot_img)

        # Create the full canvas with margins for axes
        canvas = Image.new('RGBA', (self.canvas_width, self.canvas_height), (20, 20, 25, 255))

        # Paste the plot image onto the canvas
        plot_img_rgba = plot_img.convert('RGBA')
        canvas.paste(plot_img_rgba, (self.margin_left, self.margin_top))

        # Draw axes
        canvas = self._draw_axes(canvas)

        # Convert to PNG bytes
        buffer = io.BytesIO()
        canvas.save(buffer, format='PNG')
        buffer.seek(0)

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
        ui.label('mzML Peak Map Viewer').classes('text-3xl font-bold mb-2')
        ui.label('High-performance visualization with Datashader + pyOpenMS').classes('text-gray-400 mb-4')

        # File upload section - mzML
        with ui.card().classes('w-full max-w-5xl mb-4'):
            ui.label('Load Data').classes('text-xl font-semibold mb-2')

            with ui.row().classes('w-full items-end gap-4'):
                # mzML section
                with ui.column().classes('flex-1'):
                    ui.label('mzML File (Peak Data)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        async def handle_mzml_upload(e):
                            content = e.content.read()
                            temp_path = Path('/tmp') / e.name
                            temp_path.write_bytes(content)
                            if viewer.load_mzml(str(temp_path)):
                                viewer.update_plot()

                        ui.upload(
                            label='Upload mzML',
                            on_upload=handle_mzml_upload,
                            auto_upload=True
                        ).props('accept=.mzML,.mzml').classes('w-48')

                        mzml_input = ui.input(
                            placeholder='/path/to/file.mzML'
                        ).classes('flex-1')

                        async def load_mzml_path():
                            path = mzml_input.value
                            if path and Path(path).exists():
                                if viewer.load_mzml(path):
                                    viewer.update_plot()
                            else:
                                ui.notify("File not found", type="warning")

                        ui.button('Load', on_click=load_mzml_path).props('color=primary dense')

                # FeatureXML section
                with ui.column().classes('flex-1'):
                    ui.label('FeatureXML (Feature Overlay)').classes('text-sm text-gray-400')
                    with ui.row().classes('w-full items-end gap-2'):
                        async def handle_feature_upload(e):
                            content = e.content.read()
                            temp_path = Path('/tmp') / e.name
                            temp_path.write_bytes(content)
                            if viewer.load_featuremap(str(temp_path)):
                                viewer.update_plot()

                        ui.upload(
                            label='Upload featureXML',
                            on_upload=handle_feature_upload,
                            auto_upload=True
                        ).props('accept=.featureXML,.xml').classes('w-48')

                        feature_input = ui.input(
                            placeholder='/path/to/features.featureXML'
                        ).classes('flex-1')

                        async def load_feature_path():
                            path = feature_input.value
                            if path and Path(path).exists():
                                if viewer.load_featuremap(path):
                                    viewer.update_plot()
                            else:
                                ui.notify("File not found", type="warning")

                        ui.button('Load', on_click=load_feature_path).props('color=primary dense')

                        def clear_features():
                            viewer.clear_features()
                            viewer.update_plot()

                        ui.button('Clear', on_click=clear_features).props('color=negative dense')

        # Info bar
        with ui.row().classes('w-full justify-center gap-8 mb-2'):
            viewer.info_label = ui.label('No file loaded').classes('text-gray-400')
            viewer.feature_info_label = ui.label('Features: None').classes('text-cyan-400')
            viewer.status_label = ui.label('Ready').classes('text-green-400')

        # Range display
        with ui.row().classes('w-full justify-center gap-8 mb-2'):
            viewer.rt_range_label = ui.label('RT: -- - -- s').classes('text-blue-300')
            viewer.mz_range_label = ui.label('m/z: -- - --').classes('text-blue-300')

        # Feature display options
        with ui.row().classes('w-full justify-center gap-4 mb-2'):
            ui.label('Show:').classes('text-gray-400')

            def toggle_centroids():
                viewer.show_centroids = centroid_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            centroid_cb = ui.checkbox('Centroids', value=True, on_change=toggle_centroids).classes('text-green-400')

            def toggle_bboxes():
                viewer.show_bounding_boxes = bbox_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            bbox_cb = ui.checkbox('Bounding Boxes', value=True, on_change=toggle_bboxes).classes('text-yellow-400')

            def toggle_hulls():
                viewer.show_convex_hulls = hull_cb.value
                if viewer.df is not None:
                    viewer.update_plot()

            hull_cb = ui.checkbox('Convex Hulls', value=True, on_change=toggle_hulls).classes('text-cyan-400')

        # The main plot image (now with axes included)
        with ui.card().classes('p-0'):
            viewer.image_element = ui.image().classes('w-full').style(
                f'width: {viewer.canvas_width}px; height: {viewer.canvas_height}px; background: #141419;'
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

        # Legend
        with ui.expansion('Legend & Help', icon='help').classes('w-full max-w-4xl mt-2'):
            with ui.row().classes('gap-8'):
                with ui.column():
                    ui.label('Feature Overlay Colors:').classes('font-semibold')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:#00ff64;border-radius:50%;border:1px solid white;"></div>')
                        ui.label('Centroid (feature center)')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;border:2px solid #ffff00;"></div>')
                        ui.label('Bounding Box')
                    with ui.row().classes('items-center gap-2'):
                        ui.html('<div style="width:16px;height:16px;background:rgba(0,200,255,0.5);border:1px solid #00c8ff;"></div>')
                        ui.label('Convex Hull')

                with ui.column():
                    ui.label('Keyboard Shortcuts:').classes('font-semibold')
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
