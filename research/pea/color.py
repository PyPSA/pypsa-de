colors = [
    '#002e4f', '#1b3d5b', '#2274A5', '#3f8d9e', '#44AF69',
    '#5e9a7a', '#F1C40F', '#d4a017', '#FE7F2D', '#b86f33'
]

fig, ax = plt.subplots(figsize=(10, 2))

for i, color in enumerate(colors):
    ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
    ax.text(i + 0.5, 0.5, color, ha='center', va='center', fontsize=10, color='white' if i != 5 else 'black', weight='bold')

ax.set_xlim(0, len(colors))
ax.set_ylim(0, 1)

ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])

plt.show()


import matplotlib.pyplot as plt
import matplotlib.colors as mc
import colorsys

def lighten_color(color, amount=0.2):
    try:
        c = mc.cnames[color]
    except:
        c = color
    rgb = mc.to_rgb(c)
    # Convert to HLS and adjust the L (lightness) value
    h, l, s = colorsys.rgb_to_hls(*rgb)
    l = max(0, min(1, l + amount * (1 - l)))
    return colorsys.hls_to_rgb(h, l, s)


# lower
original_colors = [
    '#002e4f', '#1b3d5b', '#2274A5', '#3f8d9e', 
    '#44AF69', '#5e9a7a', '#F1C40F', '#d4a017', 
    '#FE7F2D', '#b86f33'
]

light_colors = [lighten_color(color, 0.5) for color in original_colors]

fig, ax = plt.subplots(figsize=(10, 2))

for i, color in enumerate(light_colors):
    ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=color))
    text_color = 'black' if sum(mc.to_rgb(color)) > 1.5 else 'white'
    ax.text(i + 0.5, 0.5, original_colors[i], ha='center', va='center', fontsize=10, color=text_color, weight='bold')

ax.set_xlim(0, len(light_colors))
ax.set_ylim(0, 1)

ax.set_xticks([])
ax.set_yticks([])
ax.set_xticklabels([])
ax.set_yticklabels([])

plt.show()