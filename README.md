# NonIdealReactorApp

`NonIdealReactorApp` is a MATLAB application for teaching and studying non-ideal reactors through RTD-based analysis, reactor-performance prediction, equivalent-model interpretation, and design-oriented comparison tools. The current app covers `RTD Analysis`, direct and one-parameter prediction-model routes, and a `Design & Optimization` workspace built on top of the selected RTD, reaction system, and feed stream.

This project extends the original `ReactorApp` line developed at the University of Alicante. The current non-ideal-reactor branch keeps that educational spirit while adding RTD-driven prediction routes, visualization layers, session persistence, and a design workspace for fitting, reactive comparison, and process-side optimization.

## Requirements

- MATLAB with support for `uifigure`-based apps
- A local copy of this project extracted from the distributed ZIP package

## Installation and Launch

There is no installer. Use the project as a normal MATLAB folder:

1. Download the project ZIP.
2. Extract it to any local folder.
3. Open MATLAB.
4. In MATLAB, either:
   - set the extracted project folder as the current folder, or
   - add the extracted folder to the MATLAB path.
5. Open the `ReactorApp toolbox` folder inside MATLAB if needed.

You can launch the non-ideal app in either of these ways:

### Direct launch

Run:

```matlab
NonIdealReactorApp
```

### Launch through the original main window

Run:

```matlab
ReactorApp
```

Then open the non-ideal branch from the main window by selecting the corresponding non-ideal reactor option.

## Academic Origin and Acknowledgment

This work is rooted in the original `ReactorApp` created by **Isabela Fons** at the **University of Alicante**. The original project should be cited through:

- [Journal article](https://www.sciencedirect.com/science/article/pii/S1749772826000096?ref=pdf_download&fr=RR-2&rr=a1a6adaa1fb8c112)
- [University of Alicante repository (TFG)](https://hdl.handle.net/10045/107785)
- [MathWorks File Exchange project](https://es.mathworks.com/matlabcentral/fileexchange/76917-reactorapp-toolbox?s_tid=prof_contriblnk)

## Current NonIdealReactorApp Author

The current `NonIdealReactorApp` branch and its non-ideal-reactor extensions are authored and maintained by the present project author.

- [GitHub repository](https://github.com/JavierBerenguer/Trabajo-Final-de-Grado.git)
- MathWorks project page  [![View Trabajo-Final-de-Grado on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://es.mathworks.com/matlabcentral/fileexchange/184241-trabajo-final-de-grado)
- [University of Alicante TFG repository entry](https://hdl.handle.net/10045/169730)

## Documentation

This folder belongs to the technical documentation branch. For user-facing help inside the app, use:

- `Help > User Guide`
- `Help > Technical Guide`
