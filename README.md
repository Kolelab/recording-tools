# recording-tools
Tools for analyzing and processing current-clamp voltage recording data

Traces Processor (trcsProcessor.m) - a matlab-based tool for removing noisy traces from current-clamp recordings at a given location, resulting in an average voltage response for each current injection protocol.
To get started, prepare for import the voltage responses at a given recording location across current injections and open roundsd.m. Run.

Velometer (velometer.m) - another matlab-based tool for calculating voltage threshold as a function of the voltage waveform's temporal derivative (Kole and Stuart 2008). Run to interactively import a given action potential waveform and define the desired voltage threshold (i.e. 100 V/s). Velometer will return the time at which that threshold is breached.

intersections.m and roundsd.m: required auxiliary tools to run trcsProcessor and velometer.
