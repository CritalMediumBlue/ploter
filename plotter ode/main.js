var pointCount = 1000;
var i, r;

let prod_R3 = 0.02; // molecule sec-1 Assuming P3 is already defined
let prod_R1 = 0.15; // Assuming P1 is already defined
let deg_R = 0.2; // h-1Assuming DR is already defined
let deg_I = 0.2; // h-1Assuming DI is already defined
let deg_L = 0.6; // h-1 Assuming DL is already defined
let for_com_RI = 0.32; // #-1 h-1 or 2.73x10^4 M^-1 s^-1 Assuming KonRI is already defined
let for_com_RL = 0.32; // #-1 h-1 or 4.85 x10^3 M^-1 s^-1Assuming KonRL is already defined
let prod_L = 0.2; // Molecules second-1 Assuming PL is already defined
let n_R = 1; // Assuming nR is already defined
let n_A = 1; // Assuming nA is already defined
let K_R = 3; // Assuming KR is already defined
let K_A = 3; // Assuming KA is already defined
let A = 6;


var initialValues = [];
var grid_size = 4;
var range = [0, 6]; // Adjust as needed

for (let x = 0; x < grid_size; x++) {
    for (let y = 0; y < grid_size; y++) {
        for (let z = 0; z < grid_size; z++) {
            initialValues.push([
                range[0] + (range[1] - range[0]) * x / (grid_size - 1),
                range[0] + (range[1] - range[0]) * y / (grid_size - 1),
                range[0] + (range[1] - range[0]) * z / (grid_size - 1)
            ]);
        }
    }
}
var data = [];

for (let j = 0; j < initialValues.length; j++) {
    var R = [initialValues[j][0]];
    var I = [initialValues[j][1]];
    var L = [initialValues[j][2]];
    var c = [];

    for(i = 0; i < pointCount; i++) {
        // dR/dt equation
        let activation_R = prod_R1 / (1 + Math.pow(R[i] / K_R, n_R)) * Math.pow(A[i] / K_A, n_A) / (1 + Math.pow(A[i] / K_A, n_A));
        let dR = prod_R3 + activation_R - deg_R * R[i] - for_com_RI * R[i] * I[i] - for_com_RL * R[i] * L[i]; /* + dis_com_L * ComL[i] + dis_com_I * ComI[i]; */
        R.push(R[i] + dR);
        
        // dI/dt equation
        let dI = activation_R - deg_I * I[i] - for_com_RI * R[i] * I[i]; /* + dis_com_I * ComI[i]; */
        I.push(I[i] + dI);
        
        // dL/dt equation
        let activation_L = prod_L / (1 + Math.pow(R[i] / K_R, n_R));
        let dL = activation_L - deg_L * L[i] - for_com_RL * R[i] * L[i]; /* + dis_com_L * ComL[i]; */
        L.push(L[i] + dL);

        c.push(i)
    }

    data.push({
        type: 'scatter3d',
        mode: 'lines+markers',
        x: R, //SinR
        y: I, //SinI
        z: L, //SlrR
        line: {
            width: 5,
            color: c,
            colorscale: "Viridis",
            cmin: 0,
            cmax: 1000
        },
        marker: {
            size: 1,
            color: c,
            colorscale: "Viridis",
            cmin: 0,
            cmax: 1000
        }
    });
}

Plotly.newPlot('myDiv', data, {
    scene: {
        xaxis: {title: 'SinR'},
        yaxis: {title: 'SinI'},
        zaxis: {title: 'SlrR'}
    },
    showlegend: false
});


