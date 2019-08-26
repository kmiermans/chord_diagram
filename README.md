# chord_diagram
Makes a circos-like diagram to visualize 'chords' that connect points on the circle

![example chord diagram](https://github.com/kmiermans/chord_diagram/blob/master/example.png)

## Summary
This packages leverages matplotlib to draw a circle with chords that connect points on that circle. I designed this module to represent a circular polymer with bivalent proteins, connecting to monomers on that polymer to each other. Hence, the length of the circle (*kwarg*: polymer_length) represents the length of the polymer, and the input array (*kwarg*: sites) is an Nx2 array representing all the different connections between those points on the circle.

## Terms of use
You're free to use this script, but please cite it for use in published works as *C.A.Miermans, mayavi_tubeplot, (2018), GitHub repository, https://github.com/kmiermans/mayavi_tubeplot*.


## Usage
Usage is incredibly simple:
```
sites = np.array([[1,5], [4, 10]])

X = ProteinChordDiagram(sites, 12, color_scheme='light')
X.draw_all(add_text=True, alpha=0.5)

X.savefig('example.png')
```
and you're done!

## Installation
This class is written in Python 3. Converting to Python 2 would require some slight changes.

## Dependencies
- numpy
- matplotlib
- seaborn (for color schemes)
