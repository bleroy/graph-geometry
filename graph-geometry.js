// Graph Geometry, an experiment in emerging geometries from graphs by Bertrand Le Roy

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

const maxDimension = 6;

// Some utilities

/**
 * A simple vertex in an unoriented graph
 */
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

/**
 * A vertex that corresponds to a point in  some 3D space
 */
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

const interpolationDenominators = new Array(maxDimension);
for (let i = 0; i < maxDimension; i++) {
    let denominator = 1;
    for (let j = 0; j < maxDimension; j++) {
        if (i !== j) denominator *= i - j;
    }
    interpolationDenominators[i] = denominator;
}

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

const evaluate = (polynomial, x) => polynomial.reduce(
    (previous, coeff, i) => previous + coeff * Math.pow(x, i), 0
);

// The math stuff //

// First, prepare the graphs

// Create a bunch of random points on a sphere of radius 1
const vertexCount = 300;
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

const sphereGenerationCounts = countNeighborGenerations(vertices[0]);
const sphereLogNeighborCounts = toDimension(sphereGenerationCounts);

const sphereInterpolationPolynomial = getInterpolatedPolynomial(sphereGenerationCounts);
console.log(sphereInterpolationPolynomial);
const sphereInterpolatedValues = numbers(maxDimension + 1).map(i => evaluate(sphereInterpolationPolynomial, i));
console.log(sphereInterpolatedValues);

console.log(`Spherical random graph generation vertex counts: ${sphereGenerationCounts}`)
console.log(`Spherical random graph generation dimension estimates: ${sphereLogNeighborCounts}`)

// Simple 1D, 2D, 3D, and 4D spaces with a regular square mesh
const line = [];
const plane = [];
const volume = [];
const hypervolume = [];
const extent = 13;
for (let i = 0; i < extent; i++) {
    const newLineVertex = new Vertex();
    if (i > 0) newLineVertex.connect(line[i - 1]);
    line.push(newLineVertex);

    for (let j = 0; j < extent; j++) {
        const newPlaneVertex = new Vertex();
        if (i > 0) newPlaneVertex.connect(plane[(i - 1) * extent + j]);
        if (j > 0) newPlaneVertex.connect(plane[i * extent + j - 1]);
        plane.push(newPlaneVertex);

        for (let k = 0; k < extent; k++) {
            const newVolumeVertex = new Vertex();
            if (i > 0) newVolumeVertex.connect(volume[(i - 1) * extent * extent + j * extent + k]);
            if (j > 0) newVolumeVertex.connect(volume[i * extent * extent + (j - 1) * extent + k]);
            if (k > 0) newVolumeVertex.connect(volume[i * extent * extent + j * extent + k - 1]);
            volume.push(newVolumeVertex);

            for (let l = 0; l < extent; l++) {
                const newHyperVolumeVertex = new Vertex();
                if (i > 0) newHyperVolumeVertex.connect(hypervolume[(i - 1) * extent * extent * extent + j * extent * extent + k * extent + l]);
                if (j > 0) newHyperVolumeVertex.connect(hypervolume[i * extent * extent * extent + (j - 1) * extent * extent + k * extent + l]);
                if (k > 0) newHyperVolumeVertex.connect(hypervolume[i * extent * extent * extent + j * extent * extent + (k - 1) * extent + l]);
                if (l > 0) newHyperVolumeVertex.connect(hypervolume[i * extent * extent * extent + j * extent * extent + k * extent + l - 1]);
                hypervolume.push(newHyperVolumeVertex);
            }
        }
    }
}

const mid = Math.floor(extent / 2);
const lineGenerationCounts = countNeighborGenerations(line[mid]);
const lineLogNeighborCounts = toDimension(lineGenerationCounts);

const lineInterpolationPolynomial = getInterpolatedPolynomial(lineGenerationCounts);
console.log(lineInterpolationPolynomial);
const lineInterpolatedValues = numbers(maxDimension + 1).map(i => evaluate(lineInterpolationPolynomial, i));
console.log(lineInterpolatedValues);

console.log(`linear graph generation vertex counts: ${lineGenerationCounts}`)
console.log(`Linear graph generation dimension estimates: ${lineLogNeighborCounts}`)

const planeGenerationCounts = countNeighborGenerations(plane[mid * extent + mid]);
const planeLogNeighborCounts = toDimension(planeGenerationCounts);

const planeInterpolationPolynomial = getInterpolatedPolynomial(planeGenerationCounts);
console.log(planeInterpolationPolynomial);
const planeInterpolatedValues = numbers(maxDimension + 1).map(i => evaluate(planeInterpolationPolynomial, i));
console.log(planeInterpolatedValues);

console.log(`Planar graph generation vertex counts: ${planeGenerationCounts}`)
console.log(`Planar graph generation dimension estimates: ${planeLogNeighborCounts}`)

const volumeGenerationCounts = countNeighborGenerations(volume[mid * extent + mid]);
const volumeLogNeighborCounts = toDimension(volumeGenerationCounts);

const volumeInterpolationPolynomial = getInterpolatedPolynomial(volumeGenerationCounts);
console.log(volumeInterpolationPolynomial);
const volumeInterpolatedValues = numbers(maxDimension + 1).map(i => evaluate(volumeInterpolationPolynomial, i));
console.log(volumeInterpolatedValues);

console.log(`Volume graph generation vertex counts: ${volumeGenerationCounts}`)
console.log(`Volume graph generation dimension estimates: ${volumeLogNeighborCounts}`)

const hypervolumeGenerationCounts = countNeighborGenerations(hypervolume[mid * extent + mid]);
const hypervolumeLogNeighborCounts = toDimension(hypervolumeGenerationCounts);

const hypervolumeInterpolationPolynomial = getInterpolatedPolynomial(hypervolumeGenerationCounts);
console.log(hypervolumeInterpolationPolynomial);
const hypervolumeInterpolatedValues = numbers(maxDimension + 1).map(i => evaluate(hypervolumeInterpolationPolynomial, i));
console.log(hypervolumeInterpolatedValues);

console.log(`Hypervolume graph generation vertex counts: ${hypervolumeGenerationCounts}`)
console.log(`Hypervolume graph generation dimension estimates: ${hypervolumeLogNeighborCounts}`)

// The rendering stuff //

// Scene and camera
const scene = new Scene();
const camera = new PerspectiveCamera(75, window.innerWidth / window.innerHeight * 3 / 2, 0.1, 1000);

// Render surface
const renderer = new WebGLRenderer();
renderer.setSize( window.innerWidth, window.innerHeight * 2 / 3 );
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

// Add our random points to the scene
//const dotGeometry = new Geometry();
const lineGeometry = new Geometry();
for (let i = 0; i < vertices.length; i++) {
    const vertex = vertices[i];
//    dotGeometry.vertices.push(vertex.spatialLocation);
    for (let neighborIndex in vertex.neighbors) {
        lineGeometry.vertices.push(vertex.spatialLocation);
        lineGeometry.vertices.push(vertex.neighbors[neighborIndex].spatialLocation);
    }
}
//const dotMaterial = new PointsMaterial({ size: 1, sizeAttenuation: false });
//const dots = new Points(dotGeometry, dotMaterial);
//scene.add(dots);
const lineMaterial = new LineBasicMaterial({ color: 0x808080, linewidth: 1 });
const lines = new LineSegments(lineGeometry, lineMaterial);
scene.add(lines);

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

const chart = new Chart(graphCanvas.getContext('2d'), {
    type: 'bar',
    data: {
        labels: numbers(sphereLogNeighborCounts.length).map(i => `x^${i}`),
        datasets: [{
            label: 'Spherical graph',
            data: sphereInterpolationPolynomial,
            backgroundColor: 'blue'
        }, {
            label: '1D grid',
            data: lineInterpolationPolynomial,
            backgroundColor: 'red'
        }, {
            label: '2D grid',
            data: planeInterpolationPolynomial,
            backgroundColor: 'green'
        }, {
            label: '3D grid',
            data: volumeInterpolationPolynomial,
            backgroundColor: 'magenta'
        }, {
            label: '4D grid',
            data: hypervolumeInterpolationPolynomial,
            backgroundColor: 'cyan'
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