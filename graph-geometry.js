// Graph Geometry, an experiment in emerging geometries from graphs; (c) 2018 Bertrand Le Roy

// THREE.js imports
const {
    AmbientLight,
    Euler,
    Geometry,
    LineBasicMaterial,
    LineSegments,
    OrbitControls,
    PerspectiveCamera,
    PointLight,
    Points,
    PointsMaterial,
    Quaternion,
    Scene,
    Vector3,
    WebGLRenderer
} = THREE;

/** Setting how far we'll explore the graph when probing for dimensionality, etc. */
const maxDimension = 6;

// Some utilities

/** A simple vertex in an unoriented graph */
class Vertex {
    constructor() {
        this.neighbors = [];
    }

    /**
     * Connects to another vertex. The connection goes both ways.
     * @param {Vertex} vertex the vertex to connect to
     */
    connect(vertex) {
        if (!this.neighbors.includes(vertex)) this.neighbors.push(vertex);
        if (!vertex.neighbors.includes(this)) vertex.neighbors.push(this);
    }

    /**
     * Checks if a vertex is a direct neighbor.
     * @param {Vertex} vertex The vertex to look for
     * @returns {Boolean} true if the vertex is a neighbor
     */
    isConnectedTo(vertex) {
        return this.neighbors.includes(vertex);
    }
}

/** A vertex that corresponds to a point in  some 3D space */
class Point3 extends Vertex {
    /**
     * Contructs a Point3 from its coordinates.
     * @param {Number} x The x coordinate
     * @param {Number} y The y coordinate
     * @param {Number} z The z coordinate
     */
    constructor(x, y, z) {
        super();
        this.spatialLocation = new Vector3(x, y, z);
    }
}

/**
 * Creates an array with all integers between zero and the provided number.
 * This is useful to create more complicated arrays of numbers.
 * @example `numbers(5).map(n => n * n)` returns `[0, 1, 4, 9, 16]`.
 * @param {Number} l The number of integers to include in the array.
 * @returns {Number[]} An array containing all integers between 0 and {@link l} (not included).
 */
const numbers = l => [...Array(l).keys()];

/**
 * Returns the last item in an array, or undefined if it's empty.
 * @param {Array} array the array
 * @returns {*} the last item in array, or undefined
 */
const last = array => array.length > 0 ? array[array.length - 1] : undefined;

/**
 * Returns the square of the distance between two points.
 * @param {Vector3} p1 the first point
 * @param {Vector3} p2 the second point
 * @returns {Number} the squared distance between p1 and p2
 */
function sqdist(p1, p2) {
    const x = p1.x - p2.x, y = p1.y - p2.y, z = p1.z - p2.z;
    return x*x + y*y + z*z;
};

/** Constants used in polynomial interpolation */
const interpolationDenominators = new Array(maxDimension);
for (let i = 0; i < maxDimension; i++) {
    let denominator = 1;
    for (let j = 0; j < maxDimension; j++) {
        if (i !== j) denominator *= i - j;
    }
    interpolationDenominators[i] = denominator;
}

/** Constants used in polynomial interpolation */
const interpolationPolynomials = new Array(maxDimension);
for (let i = 0; i < maxDimension; i++) {
    interpolationPolynomials[i] = new Array(maxDimension);
    interpolationPolynomials[i].fill(0);
    interpolationPolynomials[i][0] = 1 / interpolationDenominators[i];
    for (let j = 0; j < maxDimension; j++) {
        if (j === i) continue;
        const newCoefficients = new Array(maxDimension);
        newCoefficients.fill(0);
        for (let k = 0; k < ((j < i) ? j + 1 : j); k++) {
          newCoefficients[k + 1] += interpolationPolynomials[i][k];
          newCoefficients[k] -= j * interpolationPolynomials[i][k];
       }   
       interpolationPolynomials[i] = newCoefficients;
    }
}

/**
 * Interpolates a function with a polynomial.
 * If A is the array of coefficients returned by the function, the function
 * can be evaluated using the {@link evaluate} function, passing in A and the value
 * to evaluate the function for.
 * @param {Number[]} data The values of the function to interpolate for any number
 * of consecutive integer values of x, starting with 0.
 * @return {Number[]} The coefficients of the polynomial interpolation for the provided function.
 */
function getInterpolatedPolynomial(data) {
    const polynomial = new Array(maxDimension);
    polynomial.fill(0);
    for (let i = 0; i < maxDimension; i++) {
      for (let j = 0; j < maxDimension; j++) {
         polynomial[j] += data[i] * interpolationPolynomials[i][j];
      }
    }
    return polynomial;
}

/**
 * Evaluates the value of a polynomial interpolation at a specific point.
 * @param {Number[]} polynomial The coefficients of the interpolation polynomial.
 * @param {Number} x The value to evaluate the function for.
 * @returns {Number} The interpolated value, that is the sum of `polynomial[i] * x^i`
 * for i going over all the indices in the {@link polynomial} array of coefficients.
 */
const evaluate = (polynomial, x) => polynomial.reduce(
    (previous, coeff, i) => previous + coeff * Math.pow(x, i), 0
);

/**
 * Creates and appends an HTML element.
 * @param {String} tag The element name
 * @param {Object | String} attr Key-value pairs describing the attributes, or child elements if the value is an array,
 * or a simple string if the element must be a simple tag with text contents.
 * Child elements have a single property whose name is the child tag name, and value is the set of attributes.
 * If a child element is a simple string, it is added to the content of the element.
 * If an attribute is a function, it is added as an event handler.
 * @param {HTMLElement} parent The element to append the new element to.
 * @returns the created element.
 */
function addElement(tag, attr, parent = document.body) {
    const el = document.createElement(tag);
    if (typeof(attr) === 'string') {
        const content = document.createTextNode(attr);
        el.appendChild(content);
        parent.appendChild(el);
        return el;
    }
    for (let name in attr) {
        const value = attr[name];
        if (Array.isArray(value)) {
            for (let child of value) {
                if (typeof(child) === 'string') {
                    el.appendChild(document.createTextNode(child));
                }
                else {
                    for (let childTag in child) {
                        // Should be only one of those, otherwise there be dragons.
                        addElement(childTag, child[childTag], el);
                    }
                }
            }
        }
        else {
            if (typeof(value) === 'function') {
                el.addEventListener(name, value);
            }
            else el[name] = value;
        }
    }
    parent.appendChild(el);
    return el;
}

// The math stuff //

// First, prepare the graphs

/**
 * Creates a bunch of random points on a sphere of radius 1, then connects them,
 * trying to keep connections short while connecting the whole set to create a single
 * graph.
 * @param {Number} vertexCount The number of points to randomly create.
 * @returns {Point3[]} The array of points.
 */
function createRandomSphericalMesh(vertexCount) {
    const vertices = [];
    const randomAngle = () => (Math.random() - 0.5) * 2 * Math.PI;
    let transform = new Quaternion(1, 0, 0, 0);
    for (let i = 0; i < vertexCount; i++) {
        const newDot = new Point3(1, 0, 0);
        newDot.spatialLocation.applyQuaternion(transform);
        vertices.push(newDot);
        const additionalRotation = new Euler(randomAngle(), randomAngle(), randomAngle());
        transform = transform.multiply(new Quaternion().setFromEuler(additionalRotation));
    }

    /**
     * Try to connect the whole set of points by finding the nearest spatial neighbors.
     * This is expensive, but we're just doing it once to build a test mesh.
     *
     * Start from an arbitrary point, then find the closest non-neighbor to any of the currently
     * connected set (even if it's already in the set) and add a connection.
     * Finish if all vertices have been connected.
     * @param {Point3[]} vertices the set of points to connect
     * @param {Number} startingIndex the index of the starting point in the vertices array
     * @returns {Number[]} the set of vertex indices, in the order in which they got connected
     */
    function connectClosestVertices(vertices, startingIndex = 0) {
        const connected = [startingIndex];
        while (connected.length < vertices.length) {
            let smallestSqdist = Number.POSITIVE_INFINITY;
            let closestFromIndex = -1, closestToIndex = -1;
            for (let i of connected) {
                const vertex = vertices[i];
                for (let j = 0; j < vertices.length; j++) {
                    if (i === j || vertex.isConnectedTo(vertices[j])) continue;
                    const newSqdist = sqdist(vertex.spatialLocation, vertices[j].spatialLocation);
                    //console.log(`Distance between P${i} and P${j} is ${newSqdist}`)
                    if (newSqdist < smallestSqdist) {
                        smallestSqdist = newSqdist;
                        closestFromIndex = i;
                        closestToIndex = j;
                    }
                }
            }
            if (closestToIndex === -1) break;

            const toNeighbor = vertices[closestToIndex];
            const fromNeighbor = vertices[closestFromIndex];
            fromNeighbor.connect(toNeighbor);
            //console.log(`Connecting P${closestFromIndex} to P${closestToIndex}`);
            if (!connected.includes(closestToIndex)) {
                connected.push(closestToIndex);
                //console.log(`Connected set: {${connected.map(p => `P${p}`).join(', ')}}`);
            }
        }
        return connected;
    }

    // Connect the vertices, then do it again from the last connected point, which had the least opportunity for connections
    var connected = connectClosestVertices(vertices);
    connectClosestVertices(vertices, last(connected));

    return vertices;
}

// Then forget about real geometry, and try to compute geometry from the graphs alone.

/**
 * Get all neighbors of all the vertices passed in.
 * @param {Vertex[]} vertices The set of vertices for which to find the complete set of neighbors
 * @returns {Vertex[]} the union of all neighbors
 */
function getNeighbors(vertices) {
    let result = new Set(vertices);
    for (let vertex of vertices) {
        for (let neighbor of vertex.neighbors) result.add(neighbor);
    }
    return [...result];
}

/**
 * Counts neighbors by generations, starting from direct neighbors, to an arbitrary number of generations.
 * Each generation's count includes all generations before it.
 * @param {Vertex} vertex The vertex for which to count neighbor generations
 * @param {Number} generations The number of generations to explore
 * @returns {Number[]} the number of neighbors for each generation
 */
function countNeighborGenerations(vertex, generations = maxDimension) {
    const neighborCounts = [1];
    let neighbors = [vertex];
    for (let generation = 0; generation < generations; generation++) {
        neighbors = getNeighbors(neighbors);
        neighborCounts.push(neighbors.length);
    }
    return neighborCounts;
}

/**
 * Estimates the dimension from the growth of the number of vertices per generation
 * @param {Number[]} countPerGeneration A list of vertex counts as one goes away from the initial vertex
 * @returns {Number[]} the estimate for the dimension at each generation. Count is one less than number of generations.
 */
function toDimension(countPerGeneration) {
    results = [];
    const offset = countPerGeneration[0];
    for (let i = 2; i < countPerGeneration.length; i++) {
        results.push(Math.log((countPerGeneration[i] - offset) / (countPerGeneration[i - 1] - offset)) / Math.log(i / (i - 1)));
    }
    return results;
}

const extent = maxDimension * 2 + 1;
const mid = maxDimension;

/**
 * Generates simple square meshes in 1D, 2D, 3D, and 4D.
 * @param {Number} extent The number of points along one dimension of each of the meshes.
 * @returns {Object} meshes
 * @returns {Point3[]} meshes.line The set of points on the 1D mesh.
 * @returns {Point3[]} meshes.plane The set of points on the 2D mesh.
 * @returns {Point3[]} meshes.volume The set of points on the 3D mesh.
 * @returns {Point3[]} meshes.hypervolume The set of points on the 4D mesh.
 */
function createSquareMeshes(extent) {
    const line = [];
    const plane = [];
    const volume = [];
    const hypervolume = [];
    for (let i = 0; i < extent; i++) {
        const newLineVertex = new Point3(0, 0, (i - mid) / extent);
        if (i > 0) newLineVertex.connect(line[i - 1]);
        line.push(newLineVertex);

        for (let j = 0; j < extent; j++) {
            const newPlaneVertex = new Point3(0, (j - mid) / extent, (i - mid) / extent);
            if (i > 0) newPlaneVertex.connect(plane[(i - 1) * extent + j]);
            if (j > 0) newPlaneVertex.connect(plane[i * extent + j - 1]);
            plane.push(newPlaneVertex);

            for (let k = 0; k < extent; k++) {
                const newVolumeVertex = new Point3((k - mid) / extent, (j - mid) / extent, (i - mid) / extent);
                if (i > 0) newVolumeVertex.connect(volume[(i - 1) * extent * extent + j * extent + k]);
                if (j > 0) newVolumeVertex.connect(volume[i * extent * extent + (j - 1) * extent + k]);
                if (k > 0) newVolumeVertex.connect(volume[i * extent * extent + j * extent + k - 1]);
                volume.push(newVolumeVertex);

                for (let l = 0; l < extent; l++) {
                    const displacement = l /extent;
                    const newHyperVolumeVertex = new Point3((k - mid) / extent, (j - mid) / extent, (i - mid) / extent);
                    if (i > 0) newHyperVolumeVertex.connect(hypervolume[(i - 1) * extent * extent * extent + j * extent * extent + k * extent + l]);
                    if (j > 0) newHyperVolumeVertex.connect(hypervolume[i * extent * extent * extent + (j - 1) * extent * extent + k * extent + l]);
                    if (k > 0) newHyperVolumeVertex.connect(hypervolume[i * extent * extent * extent + j * extent * extent + (k - 1) * extent + l]);
                    if (l > 0) newHyperVolumeVertex.connect(hypervolume[i * extent * extent * extent + j * extent * extent + k * extent + l - 1]);
                    hypervolume.push(newHyperVolumeVertex);
                }
            }
        }
    }
    return {line, plane, volume, hypervolume};
}

const {line, plane, volume, hypervolume} = createSquareMeshes(extent);

// The rendering stuff //

function displayVertices(vertices, startVertex, legend) {
    const generationCounts = countNeighborGenerations(startVertex);
    
    const interpolationPolynomial = getInterpolatedPolynomial(generationCounts);
    // console.log(interpolationPolynomial);
    // const interpolatedValues = numbers(maxDimension + 1).map(i => evaluate(interpolationPolynomial, i));
    // console.log(interpolatedValues);
    
    // console.log(`Graph generation vertex counts: ${generationCounts}`)

    // const logNeighborCounts = toDimension(generationCounts);
    // console.log(`Graph generation dimension estimates: ${logNeighborCounts}`)
    drawGraph(vertices);
    chart(interpolationPolynomial, legend);
}

function displayRandomSpherical() {
    const vertices = createRandomSphericalMesh(300);
    const vertex = vertices[Math.floor(Math.random() * vertices.length)];
    displayVertices(vertices, vertex, "Random Spherical");
}

function display1D() {
    displayVertices(line, line[mid], "1D mesh");
}

function display2D() {
    displayVertices(plane, plane[mid * extent + mid], "2D mesh");
}

function display3D() {
    displayVertices(volume, volume[mid * extent * extent + mid * extent + mid], "3D mesh");
}

function display4D() {
    console.log("4D");
}

// Toolbar

const toolbarHeight = 100;

const toolbar = addElement("fieldset", {
    width: window.innerWidth,
    height: toolbarHeight,
    _: [
        {legend: "Pick a graph type:"},
        {button: {
            _: ["Random spherical"],
            click: displayRandomSpherical
        }},
        " ",
        {button: {
            _: ["1D mesh"],
            click: display1D
        }},
        " ",
        {button: {
            _: ["2D mesh"],
            click: display2D
        }},
        " ",
        {button: {
            _: ["3D mesh"],
            click: display3D
        }},
        " ",
        {button: {
            _: ["4D mesh"],
            click: display4D
        }}
    ]
});

// Scene and camera
const scene = new Scene();
const camera = new PerspectiveCamera(75, window.innerWidth / (window.innerHeight - toolbarHeight) * 3 / 2, 0.1, 1000);

// Render surface
const renderer = new WebGLRenderer();
renderer.setSize( window.innerWidth, window.innerHeight * 2 / 3 - toolbarHeight );
document.body.appendChild( renderer.domElement );

// Lights
const lights = [new PointLight(0xffffff, 1, 0), new PointLight(0xffffff, 1, 0), new PointLight(0xffffff, 1, 0)];

lights[0].position.set(0, 200, 0);
lights[1].position.set(100, 200, 100);
lights[2].position.set(-100, -200, -100);

scene.add(new AmbientLight(0x404040));
scene.add(lights[0]);
scene.add(lights[1]);
scene.add(lights[2]);

function drawGraph(vertices) {
    function edgeId(id1, id2) {
        return id1 < id2 ? id1 + "-" + id2 : id2 + "-" + id1;
    }
    const oldLines = scene.getObjectByName("lines");
    if (oldLines) scene.remove(oldLines);
    // Add our random points to the scene
    //const dotGeometry = new Geometry();
    const addedEdges = {};
    const lineGeometry = new Geometry();
    for (let i = 0; i < vertices.length; i++) {
        const vertex = vertices[i];
    //    dotGeometry.vertices.push(vertex.spatialLocation);
        for (let neighborIndex in vertex.neighbors) {
            const id = edgeId(i, neighborIndex);
            if (!addedEdges[id]) {
                lineGeometry.vertices.push(vertex.spatialLocation);
                lineGeometry.vertices.push(vertex.neighbors[neighborIndex].spatialLocation);
                addedEdges[id] = true;
            }
        }
    }
    //const dotMaterial = new PointsMaterial({ size: 1, sizeAttenuation: false });
    //const dots = new Points(dotGeometry, dotMaterial);
    //scene.add(dots);
    const lineMaterial = new LineBasicMaterial({ color: 0x808080, linewidth: 1 });
    const lines = new LineSegments(lineGeometry, lineMaterial);
    lines.name = "lines";
    scene.add(lines);
}

// Camera handling, controls, and animation
camera.position.z = 2;

const controls = new OrbitControls(camera);
controls.update();

const animate = function () {
    requestAnimationFrame(animate);

    controls.update();

    renderer.render(scene, camera);
};

animate();

// Graph results
const graphCanvas = document.createElement("canvas");
graphCanvas.width = window.innerWidth;
graphCanvas.height = window.innerHeight / 3;
document.body.appendChild(graphCanvas);

function chart(data, legend) {
    new Chart(graphCanvas.getContext('2d'), {
        type: 'bar',
        data: {
            labels: numbers(data.length).map(i => `x^${i}`),
            datasets: [{
                label: legend,
                data,
                backgroundColor: 'purple'
            }]
        },
        options: {
            responsive: true,
            scales: {
                xAxes: [{
                    display: true,
                    scaleLabel: {
                        display: true,
                        labelString: 'Power of generation'
                    }
                }],
                yAxes: [{
                    display: true,
                    scaleLabel: {
                        display: true,
                        labelString: 'Coefficient'
                    }
                }]
            },
            title: {
                display: true,
                text: 'Polynomial interpolation of the number of neighbors by generation'
            }
        }
    });
}