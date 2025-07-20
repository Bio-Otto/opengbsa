#!/usr/bin/env python3
"""
OpenGBSA Logo Generator
Creates a flat, minimalist, tech-science logo for the OpenGBSA package
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, Rectangle, Polygon, RegularPolygon
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import os

def create_hexagonal_grid(ax, x_min, x_max, y_min, y_max, hex_size=0.3, alpha=0.1):
    """Create a faint hexagonal grid reminiscent of graphene"""
    hex_radius = hex_size
    hex_height = hex_radius * np.sqrt(3)
    
    for x in np.arange(x_min, x_max, hex_radius * 1.5):
        for y in np.arange(y_min, y_max, hex_height):
            # Offset every other row
            x_offset = hex_radius * 0.75 if int((y - y_min) / hex_height) % 2 == 1 else 0
            hex_x = x + x_offset
            
            # Create hexagon
            hexagon = RegularPolygon((hex_x, y), 6, radius=hex_radius, 
                                   facecolor='none', edgecolor='#333333', 
                                   linewidth=0.5, alpha=alpha)
            ax.add_patch(hexagon)

def create_benzene_ring(ax, x, y, size=0.4, color='#00E5FF'):
    """Create a hexagonal benzene ring"""
    # Main hexagon
    hexagon = RegularPolygon((x, y), 6, radius=size, 
                           facecolor='none', edgecolor=color, 
                           linewidth=2, alpha=0.8)
    ax.add_patch(hexagon)
    
    # Inner circle representing aromatic ring
    circle = Circle((x, y), size * 0.6, facecolor='none', 
                   edgecolor=color, linewidth=1, alpha=0.6)
    ax.add_patch(circle)
    
    # Add small circles at vertices (carbon atoms)
    for i in range(6):
        angle = i * np.pi / 3
        atom_x = x + size * np.cos(angle)
        atom_y = y + size * np.sin(angle)
        atom = Circle((atom_x, atom_y), size * 0.1, 
                     facecolor=color, edgecolor='none', alpha=0.8)
        ax.add_patch(atom)

def create_water_drop(ax, x, y, size=0.3, color='#00E5FF'):
    """Create a smooth water-drop icon"""
    # Water drop shape using polygon
    drop_points = [
        (x, y + size * 0.8),  # Top
        (x - size * 0.3, y + size * 0.4),  # Left curve
        (x - size * 0.2, y - size * 0.2),  # Left bottom
        (x, y - size * 0.3),  # Bottom
        (x + size * 0.2, y - size * 0.2),  # Right bottom
        (x + size * 0.3, y + size * 0.4),  # Right curve
    ]
    
    drop = Polygon(drop_points, facecolor=color, edgecolor='none', alpha=0.6)
    ax.add_patch(drop)

def create_starburst(ax, x, y, size=0.15, color='#00E5FF'):
    """Create a tiny starburst accent"""
    # Create star shape
    star_points = []
    for i in range(10):
        angle = i * np.pi / 5
        radius = size if i % 2 == 0 else size * 0.5
        star_x = x + radius * np.cos(angle)
        star_y = y + radius * np.sin(angle)
        star_points.append((star_x, star_y))
    
    star = Polygon(star_points, facecolor=color, edgecolor='none', alpha=0.8)
    ax.add_patch(star)

def create_svg_logo():
    """Create flat, minimalist SVG logo for OpenGBSA"""
    
    # Create figure with 1:1 square aspect ratio
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(0, 8)
    ax.set_ylim(0, 8)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Deep matte black background
    rect = Rectangle((0, 0), 8, 8, facecolor='#0A0A0A', edgecolor='none')
    ax.add_patch(rect)
    
    # Faint hexagonal grid
    create_hexagonal_grid(ax, 0, 8, 0, 8, hex_size=0.4, alpha=0.05)
    
    # Centered wordmark "OpenGBSA"
    # Create gradient effect by drawing multiple layers
    text = "OpenGBSA"
    center_x, center_y = 4, 4.5
    
    # Electric-cyan to white gradient effect
    for i, alpha in enumerate(np.linspace(0.3, 1.0, 5)):
        color_intensity = int(255 * alpha)
        color = f'#{color_intensity:02x}{color_intensity:02x}{255:02x}'
        
        ax.text(center_x, center_y + i * 0.02, text, fontsize=24, color=color, 
                weight='bold', ha='center', family='sans-serif', alpha=0.8)
    
    # Main text
    ax.text(center_x, center_y, text, fontsize=24, color='#00E5FF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Benzene ring integrated with letter G (positioned where G would be)
    benzene_x = center_x - 1.2  # Approximate position of G
    benzene_y = center_y - 0.1
    create_benzene_ring(ax, benzene_x, benzene_y, size=0.4, color='#00E5FF')
    
    # Water drop integrated with benzene ring
    water_x = benzene_x + 0.2
    water_y = benzene_y - 0.3
    create_water_drop(ax, water_x, water_y, size=0.25, color='#00E5FF')
    
    # Starburst accent on upper-right of benzene ring
    starburst_x = benzene_x + 0.5
    starburst_y = benzene_y + 0.4
    create_starburst(ax, starburst_x, starburst_y, size=0.12, color='#00E5FF')
    
    # Subtitle
    ax.text(center_x, center_y - 1.2, 'MM/GBSA Analysis', fontsize=10, 
            color='#666666', ha='center', family='sans-serif')
    
    # Save as SVG
    plt.savefig('logo.svg', format='svg', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none')
    plt.close()
    
    print("Flat minimalist SVG logo created: logo.svg")

def create_png_logo():
    """Create flat, minimalist PNG logo for OpenGBSA"""
    
    # Create figure with 1:1 square aspect ratio
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Deep matte black background
    rect = Rectangle((0, 0), 10, 10, facecolor='#0A0A0A', edgecolor='none')
    ax.add_patch(rect)
    
    # Faint hexagonal grid
    create_hexagonal_grid(ax, 0, 10, 0, 10, hex_size=0.5, alpha=0.05)
    
    # Centered wordmark "OpenGBSA"
    text = "OpenGBSA"
    center_x, center_y = 5, 5.5
    
    # Electric-cyan to white gradient effect
    for i, alpha in enumerate(np.linspace(0.3, 1.0, 5)):
        color_intensity = int(255 * alpha)
        color = f'#{color_intensity:02x}{color_intensity:02x}{255:02x}'
        
        ax.text(center_x, center_y + i * 0.025, text, fontsize=30, color=color, 
                weight='bold', ha='center', family='sans-serif', alpha=0.8)
    
    # Main text
    ax.text(center_x, center_y, text, fontsize=30, color='#00E5FF', 
            weight='bold', ha='center', family='sans-serif')
    
    # Benzene ring integrated with letter G
    benzene_x = center_x - 1.5
    benzene_y = center_y - 0.1
    create_benzene_ring(ax, benzene_x, benzene_y, size=0.5, color='#00E5FF')
    
    # Water drop integrated with benzene ring
    water_x = benzene_x + 0.25
    water_y = benzene_y - 0.4
    create_water_drop(ax, water_x, water_y, size=0.3, color='#00E5FF')
    
    # Starburst accent
    starburst_x = benzene_x + 0.6
    starburst_y = benzene_y + 0.5
    create_starburst(ax, starburst_x, starburst_y, size=0.15, color='#00E5FF')
    
    # Subtitle
    ax.text(center_x, center_y - 1.5, 'MM/GBSA Analysis', fontsize=12, 
            color='#666666', ha='center', family='sans-serif')
    
    # Save as PNG
    plt.savefig('logo.png', format='png', dpi=300, bbox_inches='tight', 
                facecolor='none', edgecolor='none', transparent=True)
    plt.close()
    
    print("Flat minimalist PNG logo created: logo.png")

def create_simple_logo():
    """Create a flat minimalist text-based logo"""
    
    # Create image with 1:1 square aspect ratio
    width, height = 800, 800
    img = Image.new('RGBA', (width, height), (10, 10, 10, 255))  # Deep matte black
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 48)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    except:
        font_large = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw flat minimalist OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 50
    
    # Electric-cyan color
    draw.text((x, y), text, font=font_large, fill=(0, 229, 255, 255))
    
    # Subtitle
    subtitle = "MM/GBSA Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_small)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 20
    
    draw.text((x_sub, y_sub), subtitle, font=font_small, fill=(102, 102, 102, 255))
    
    # Save
    img.save('logo_simple.png', 'PNG')
    print("Flat minimalist simple logo created: logo_simple.png")

def create_github_banner():
    """Create a flat minimalist GitHub banner"""
    
    # Create image with GitHub banner dimensions
    width, height = 1280, 640
    img = Image.new('RGBA', (width, height), (10, 10, 10, 255))  # Deep matte black
    draw = ImageDraw.Draw(img)
    
    try:
        font_large = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 72)
        font_medium = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 32)
        font_small = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 24)
    except:
        font_large = ImageFont.load_default()
        font_medium = ImageFont.load_default()
        font_small = ImageFont.load_default()
    
    # Draw flat minimalist OpenGBSA text
    text = "OpenGBSA"
    bbox = draw.textbbox((0, 0), text, font=font_large)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]
    
    x = (width - text_width) // 2
    y = (height - text_height) // 2 - 60
    
    # Electric-cyan color
    draw.text((x, y), text, font=font_large, fill=(0, 229, 255, 255))
    
    # Subtitle
    subtitle = "Molecular Mechanics / Generalized Born Surface Area Analysis"
    bbox_sub = draw.textbbox((0, 0), subtitle, font=font_medium)
    sub_width = bbox_sub[2] - bbox_sub[0]
    
    x_sub = (width - sub_width) // 2
    y_sub = y + text_height + 20
    
    draw.text((x_sub, y_sub), subtitle, font=font_medium, fill=(102, 102, 102, 255))
    
    # Feature list
    features = [
        "• MM/GBSA Energy Calculation",
        "• Per-Residue Decomposition", 
        "• Protein-Ligand Analysis",
        "• Advanced Reporting"
    ]
    
    y_features = y_sub + 80
    for i, feature in enumerate(features):
        draw.text((width//2 - 200, y_features + i * 40), feature, 
                 font=font_small, fill=(102, 102, 102, 255))
    
    # Save
    img.save('github_banner.png', 'PNG')
    print("Flat minimalist GitHub banner created: github_banner.png")

def main():
    """Generate all flat minimalist logo variants"""
    print("Generating flat minimalist OpenGBSA logos...")
    
    # Create different logo variants
    create_svg_logo()
    create_png_logo()
    create_simple_logo()
    create_github_banner()
    
    print("\nFlat minimalist logo generation complete!")
    print("Generated files:")
    print("  • logo.svg - Flat minimalist vector logo")
    print("  • logo.png - Flat minimalist high-resolution logo")
    print("  • logo_simple.png - Flat minimalist simple text logo")
    print("  • github_banner.png - Flat minimalist GitHub repository banner")
    
    print("\nUsage:")
    print("  • Add logo.png to your README.md")
    print("  • Use github_banner.png as your repository banner")
    print("  • Use logo.svg for web applications")

if __name__ == "__main__":
    main() 