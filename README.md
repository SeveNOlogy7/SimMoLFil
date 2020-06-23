# SimMoLFil
A **Mode-Locked Fiber Laser (MoLFil)** simulator.

## General info
This project simulates pulses generation and evolution in mode-locke fiber lasers.
The pulse evolution in fiber is calcuated by solving GNLSE using Interaction Picture Method Algorithm [[1](https://ieeexplore.ieee.org/abstract/document/6153336)].
The codes in folder SimMLFL were originally written for this paper [[2](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-25-4-4414)].
This project currently contains only a minimum amount of the whole codes used in [[2](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-25-4-4414)].

The foler SimMLFL contains codes described above. The foler SimMolFil contains works in progress.

## Technologies
The project is created with:
* MATLAB version: 2016b

## Examples
mainP1.m: run a simulation and produce results of Point 1 of Figure 2 in [[2](https://www.osapublishing.org/oe/abstract.cfm?uri=oe-25-4-4414)].

## Works in progress
A simple data flow graph (DFG) presentation of Mode-Locked Fiber Laser Models to seperate
layout and connection of optical components from their implementation.

A DFG presentation can be built into some evaluatable object by a Builder.

Imaging you can simply write somthing like
```MATLAB
model = defind some model use DFG
evaluaion = model.build(Some Builder)
evaluaion(Some simulation settings)
```

In this way, it is easy to evaluate a same model with 
* different simulation settings (more points, fine step, etc.),
* different algorithms to simulate some specific kind of optical components (better gain model, multimode pulse simulation, etc.),
* and even tools written in other language.

All these can be done by adding a new Builder class to generate the evaluatable.

Also, it turns out to be clear that, the idea of DFG presentation is not only suitable for Mode-locked fiber laser models, but also for Cascaded Fiber Supercontinuum Generation models.

### Examples of current DFG design
Following is codes that generate a Cascaded Fiber Supercontinuum model.

```MATLAB
% Define a fiber component
fiber = component.Fiber();

% Create some "Node"s that performs Fiber operation
% SMF stands for single mode fiber
SMF1 = simulation.Fiber("SMF1", fiber);
SMF2 = simulation.Fiber("SMF2", fiber);
SMF3 = simulation.Fiber("SMF3", fiber);

% Create a input Node
In = simulation.Input("In");

% Now simply put the Nodes together
% Think it like:
% Out = In -> SMF1 -> SMF2 -> SMF3
% It is injecting a input (like a powerful laser pulse) through three cascaded fiber sections
% No "->" operator in Matlab, so it end up like this
Out = In + SMF1 + SMF2 + SMF3;

% Generate some human-readable statement presentation of the model
statements = Out.statements()
```

The statements will contain a multiline string:
```MATLAB
t0 = In.Input()
t1 = SMF1.Fiber[fiber=component.Fiber](t0)
t2 = SMF2.Fiber[fiber=component.Fiber](t1)
t3 = SMF3.Fiber[fiber=component.Fiber](t2)
```

These staments show what any builders of the DFG model should implement to make the DFG an evaluable.

### Codes Arrangement
* **simulation package** contains codes that help reuse a few MATLAB syntax to present a MoLFil Model in DFG
and simply build the MoLFil Model to 
    * Human readable statements
    * Evaluation object, given a specific Builder implementation

* **component package** contains optical components packed in classes. 
And some helper methods to calculate parameters for those components.

### Todos
- [x] DFG generation
- [x] DFG build process
- [ ] DFG Evaluation Class (working on it)
- [ ] A pure MATLAB Builder implementation using existing codes in SimMLFL (working on it)
- [ ] Fill component package with concrete classes
- [ ] Examples
