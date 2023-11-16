var pointCount = 5000;
var i, r;

var x = [];
var y = [];
var z = [];
var c = [];
var layout = {
    width: 1500,
    height: 1500,
    scene: {
        xaxis: {title: 'Red', range: [0, 1], autorange: false},
        yaxis: {title: 'Green', range: [0, 1], autorange: false},
        zaxis: {title: 'Blue', range: [0, 1], autorange: false}
    }
};

var dx = 0.5, dy = 0.5, dz = 0.5;
for(i = 0; i < pointCount; i++) 
{
    dx = Math.abs((dx + (Math.random() - 0.5) * 0.05) % 1);
    dy = Math.abs((dy + (Math.random() - 0.5) * 0.05) % 1);
    dz = Math.abs((dz + (Math.random() - 0.5) * 0.05) % 1);

    if(dx < 0.02 || dx > 0.98 || dy < 0.02 || dy > 0.98 || dz < 0.02 || dz> 0.98) {
        break;
    }

    x.push(dx);
    y.push(dy);
    z.push(dz);
    c.push('rgb(' + Math.round(dx * 255) + ',' + Math.round(dy * 255) + ',' + Math.round(dz * 255) + ')');
}

Plotly.newPlot('myDiv', [{
    type: 'scatter3d',
    mode: 'lines+markers',
    x: x,
    y: y,
    z: z,
    line: {
        width: 5,
        color: c,
        colorscale: undefined
    },
    marker: {
        size: 3.5,
        color: c,
        colorscale: undefined
    }
}], {
    scene: {
        xaxis: {title: 'Red', range: [0, 1], autorange: false},
        yaxis: {title: 'Green', range: [0, 1], autorange: false},
        zaxis: {title: 'Blue', range: [0, 1], autorange: false}
    }
}, layout);