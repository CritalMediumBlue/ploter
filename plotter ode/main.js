var pointCount = 1000;
var i, r;

const prod_R= 0.01;  // production rate of SinR, units: nM/min or similar
const prod_I= 0.01; // production rate of SinI, units: nM/min or similar
const prod_L= 0.01; // production rate of SlrR, units: nM/min or similar
const deg_R= 0.001; // degradation rate of SinR, units: min^-1
const deg_I= 0.001; // degradation rate of SinI, units: min^-1
const deg_L= 0.001; // degradation rate of SlrR, units: min^-1
const for_com_I= 0.001; // rate of complex formation between SinR and SinI, units: nM^-1 min^-1
const for_com_L= 0.001; // rate of complex formation between SinR and SlrR, units: nM^-1 min^-1
const dis_com_I= 0.0001; // rate of complex dissociation between SinR and SinI, units: min^-1
const dis_com_L= 0.0001;  // rate of complex dissociation between SinR and SlrR, units: min^-1
const deg_com_I= 0.001; // degradation rate of complex between SinR and SinI, units: min^-1
const deg_com_L= 0.001; // degradation rate of complex between SinR and SlrR, units: min^-1
const k= 10; // affinity or threshold constant, units might vary
const n= 2; // Hill coefficient, dimensionless


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
    var ComI = [0];
    var ComL = [0];
    var c = [];

    for(i = 0; i < pointCount; i++) {
        let dR=   prod_R - deg_R*R[i]- for_com_I*R[i]*I[i]- for_com_L*R[i]*L[i] + dis_com_L*ComL[i]+ dis_com_I*ComI[i];
        R.push( R[i] +dR); //SinR

        let dI= prod_I*(k^n)/((k^n)+(R[i]^n)) - deg_I*I[i] - for_com_I*R[i]*I[i]+ dis_com_I*ComI[i] ;
        I.push( I[i] +dI);//SinI

        let dL= prod_L*(k^n)/((k^n)+(R[i]^n)) - deg_L*L[i]- for_com_L*R[i]*L[i] + dis_com_L*ComL[i];
        L.push(L[i] +dL); //SlrR

        let dComI= for_com_I*R[i]*I[i] - dis_com_I*ComI[i]- deg_com_I*ComI[i]; 
        ComI.push(ComI[i]+dComI); //Complex between SinR and SinI

        let dComL= for_com_L*R[i]*L[i] - deg_com_L*ComL[i] - dis_com_L*ComL[i];
        ComL.push(ComL[i] +dComL); //Complex between SinR and SlrR

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


