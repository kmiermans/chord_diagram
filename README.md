# chord_diagram
Makes a circos-like diagram to visualize 'chords' that connect points on the circle

![example chord diagram](https://github.com/kmiermans/chord_diagram/blob/master/example.png)

## Summary

## Usage
Usage is incredibly simple:
```
sites = [[1,5], [4, 10]]
X = TubePlot(sites, polymer_length=12, color_scheme='light')
X.draw_all()
X.savefig('my_first_chord_diagram.pdf')
```
and you're done!

## Installation
This class is written in Python 3. Converting to Python 2 would require some slight changes.

## Dependencies
- numpy
- matplotlib
- seaborn (for color schemes)
