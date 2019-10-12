# buildGrapheneFlake

Build graphene monolayer flakes of a specified radius. The algorithm is quite slow for larger flakes, but should work fine for most use cases.

## How to Use
This program uses command line options to specify the size (in terms of ring units) of the graphene flake and name of output XYZ:

```bash
python buildGrapheneFlake -r radius -f fileName
```
