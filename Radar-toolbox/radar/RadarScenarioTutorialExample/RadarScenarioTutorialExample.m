%% Radar Scenario Tutorial
% This example shows how to construct and visualize a simple radar scenario
% programmatically using the |radarScenario|, |theaterPlot|, and
% |radarDataGenerator| objects.

% Copyright 2020 The MathWorks, Inc.

%% Scenario Setup
% To begin, create an empty radar scenario.  All scenario properties have
% default values.  The scenario does not contain any platform by default.

scenario = radarScenario

%% Adding Platforms
% A scenario is comprised of objects, called platforms, upon which you may
% mount sensors and emitters.  To add a platform, use the |platform| object
% function.  Here you create a simple tower.

tower = platform(scenario)

%% Platform Identification
% When you first construct a platform, it has a |PlatformID|, which is a
% unique identifier you can use to identify the platform.  The scenario
% assigns platform identifiers in the order that the platforms are created.
% You can specify a |ClassID| to denote the platform's classification.  
% For example, here you use a 3 to denote a tower.

tower.ClassID = 3;

%% Platform Signatures
% Sensors can detect platforms.  For radar sensors, a relevant signature is
% the radar cross-section (RCS).  By default a uniform RCS signature is
% used:
tower.Signatures{1}

%%
% Load predefined cylinder RCS data and use the data to define the RCS of
% the tower platform.

load('RCSCylinderExampleData.mat','cylinder');

tower.Signatures{1} = rcsSignature('Pattern', cylinder.RCSdBsm, ...
        'Azimuth', cylinder.Azimuth, 'Elevation', cylinder.Elevation, ...
        'Frequency', cylinder.Frequency);

%% Platform Dimensions
% By default, platforms have no dimensions and are modeled as point
% targets.  You may optionally specify the length, width, and height to
% denote the extent of the object, along with an offset of the body frame
% origin from its from its geometric center. You can specify platform
% dimensions using the |Dimensions| property.
%
% You specify the dimensions of the tower like this:

tower.Dimensions = struct( ...
    'Length', 10, ...
    'Width', 10, ...
    'Height', 60, ...
    'OriginOffset', [0 0 30]);

%%
% The tower has a length, width, and height of 10, 10, and 60 meters.  The
% origin offset of [0 0 30] indicates that its body frame origin
% (rotational center) is 30 meters in the positive z-direction of the
% platform's local body axis.

%% Platform Trajectories
% You can obtain a platform's current position and orientation through its
% |Position| and |Orientation| properties.  You can obtain more information
% about platforms using the scenario's platformPoses method:
 
scenario.platformPoses

%%
% You can specify a platform's position and orientation over time using
% its |Trajectory| property.  You can specify the trajectory of a platform
% using a |kinematicTrajectory|, |waypoinTrajectory|, or |geoTrajectory|
% object.
%
% By default, a platform consists of a static |kineticTrajectory| that
% whose body axis is perfectly centered and aligned with the scenario axes:

tower.Trajectory

%%
% To obtain a pitch angle of 4 degrees, you use the |Orientation| property
% of the trajectory object.  Specify the orientation using a quaternion
% obtained from Euler angles.

tYaw = 0;
tPitch = 4;
tRoll = 0;
tower.Trajectory.Orientation = quaternion([tYaw tPitch tRoll],'eulerd','zyx','frame');

%% Axes Conventions
% Most examples in Radar Toolbox,  use a "North-East-Down" convention.
% This means that the x-axis points towards north, the y-axis points toward
% east, and the z-axis points downwards.  Also, the x-, y-, and z-
% directions of the local body frame are the forward, right, and downward
% directions, respectively.

%% Visualizing a Scenario
% The |theaterPlot| object provides an interface to plot objects
% dynamically in a 3-D scene. You may use standard MATLAB axes plotting
% methods to add or remove annotations to the theater plot's axes, which
% you can obtain via its |Parent| property.
%
% Use a |platformPlotter| to plot platforms.
%
% As expected, the tower is centered at the origin in NED coordinates with
% a pitch of 4 degrees.

tPlot = theaterPlot('XLim',[-50 50],'YLim',[-50 50],'ZLim',[-100 0]);
pPlotter = platformPlotter(tPlot,'DisplayName','tower');

pose = platformPoses(scenario);
towerPose = pose(1);
towerDims = tower.Dimensions;

plotPlatform(pPlotter, towerPose.Position, towerDims, towerPose.Orientation);
set(tPlot.Parent,'YDir','reverse', 'ZDir','reverse');
view(tPlot.Parent, -37.5, 30)

%% Adding Sensors to Platforms
% To add a radar sensor to the platform, you can add a |radarDataGenerator|
% object on the platform by using its |Sensors| property.
%
% Keep in mind that in a NED coordinate system, the "Z" direction points
% down.  Therefore, if you want to mount a radar at the top of the tower,
% you should set its "z" position to -60 meters.
%
% The |radarDataGenerator| has the option to report detections in scenario
% coordinates.  Reporting detections in scenario coordinates makes it
% easier to compare the generated detections with the positions of the
% objects in the scenario.

towerRadar = radarDataGenerator('SensorIndex', 1, ...
    'UpdateRate', 10, ...
    'MountingLocation', [0 0 -60], ...
    'ScanMode', 'No scanning', ...
    'HasINS', true, ...
    'TargetReportFormat', 'Detections', ...
    'DetectionCoordinates', 'Scenario');

tower.Sensors = towerRadar;


%% Visualizing Coverage Areas
% To see sensor coverages in a scenario, you use a coverage plotter and
% plot the coverage configuration of the scenario.  You can widen the
% theater plot's range by adjusting the limits of its internal axes:

tPlot.XLimits = [-5000 5000];
tPlot.YLimits = [-5000 5000];
tPlot.ZLimits = [-1000 0];

covPlotter = coveragePlotter(tPlot,'DisplayName','Sensor Coverage');
plotCoverage(covPlotter, coverageConfig(scenario));

%% Platform Signatures
% You can add other platforms in the scenario and adjust parameters that
% affect how other sensors observe the platforms.  You can use an
% |rcsSignature| to model what the radar mounted on the tower would see.  
%
% The following code creates a helicopter and sets its radar cross section
% omnidirectionally to a value of 40 dBsm.

helicopter = platform(scenario);
helicopter.Dimensions = struct( ...
    'Length',30, ...
    'Width', .1, ...
    'Height', 7, ...
    'OriginOffset',[0 8 -3.2]);

helicopter.Signatures = rcsSignature( ...
    'Pattern',[40 40; 40 40], ...
    'Elevation',[-90; 90], ...
    'Azimuth',[-180,180]);

%%
% You can mount more than one sensor to any platform by placing the sensors
% in a cell array before assigning to the |Sensors| property.  

helicopter.Sensors = { ...
    radarDataGenerator('SensorIndex', 2, ...
                       'UpdateRate', 20, ...
                       'MountingLocation', [22 0 0], ...
                       'MountingAngles', [-5 0 0], ...
                       'ScanMode', 'No scanning', ...
                       'HasINS', true, ...
                       'TargetReportFormat', 'Detections', ...
                       'DetectionCoordinates', 'Scenario'), ...
    radarDataGenerator('SensorIndex', 3, ...
                       'UpdateRate', 30, ...
                       'MountingLocation', [22 0 0], ...
                       'MountingAngles', [5 0 0], ...
                       'ScanMode', 'No scanning', ...
                       'HasINS', true, ...
                       'TargetReportFormat', 'Detections', ...
                       'DetectionCoordinates', 'Scenario')};
                       
%% Platform Motion and Animation
% You can arrange for the helicopter to cross the path of the radar beam.
% This shows how to make the helicopter follow a straight 100-meter path at
% a constant velocity with an elevation of 250 meters for seven seconds:

helicopter.Trajectory = waypointTrajectory([2000 50 -250; 2000 -50 -250],[0 7]);


%%
% Platform motion across time is performed by using a while-loop and
% calling the scenario's |advance| method.
%
% You can plot all the platforms positions, orientations and dimensions
% in the loop:

profiles = platformProfiles(scenario);
dimensions = vertcat(profiles.Dimensions);

while advance(scenario)
    poses = platformPoses(scenario);
    positions = vertcat(poses.Position);
    orientations = vertcat(poses.Orientation);
    plotPlatform(pPlotter, positions, dimensions, orientations);
    plotCoverage(covPlotter, coverageConfig(scenario));
    % to animate more slowly uncomment the following line
    % pause(0.01)
end

%% Detecting platforms
% In the example above, you added three radars with different update
% rates:  the tower has an update rate of 10 Hz, and the helicopter had two
% radars with update rates of 20 Hz and 30 Hz, respectively.  
%
% The scenario can be placed into a mode in which the call to |advance|
% updates the time of simulation as needed to update each of the sensors it
% contains.  You can achieve this by setting the |UpdateRate| of the
% scenario to zero.

scenario.UpdateRate = 0;

%%
% To show the simulation time, add a UI control to the figure.

fig = ancestor(tPlot.Parent,'Figure');
timeReadout = uicontrol(fig,'Style','text','HorizontalAlignment','left','Position',[0 0 200 20]);
timeReadout.String = "SimulationTime: " + scenario.SimulationTime;

%%
% Now the proper sensor times can be reached.  You can use |detect| to 
% get the detections available by each sensor within the loop.  Detections
% can be shown by constructing a |detectionPlotter| object.

dPlotter = detectionPlotter(tPlot,'DisplayName','detections');

%%
% You can run the same simulation again by restarting it and modifying it
% to report detections.
%
% The detected positions for a |radarDataGenerator| can be extracted from
% the detection |Measurement| field:

restart(scenario);

while advance(scenario)
    timeReadout.String = "SimulationTime: " + scenario.SimulationTime;
    detections = detect(scenario);

    % extract column vector of measurement positions
    allDetections = [detections{:}];
    if ~isempty(allDetections)
        measurement = cat(2,allDetections.Measurement)';
    else
        measurement = zeros(0,3);
    end

    plotDetection(dPlotter, measurement);
    
    poses = platformPoses(scenario);
    positions = vertcat(poses.Position);
    orientations = vertcat(poses.Orientation);
    plotPlatform(pPlotter, positions, dimensions, orientations);
    plotCoverage(covPlotter, coverageConfig(scenario));
end

%%
% Notice that the update time of simulation increments non-uniformly to the
% times required by each of the sensors.

%% Summary
% In this example, you learned how to construct and visualize a simple
% scenario and obtain detections generated by a radar data generator.