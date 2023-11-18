var pointCount = 50000; //5000 = 2.7 hours = 166.7 minutes
var i, r;

/* - Weakly transcribed genes: 0.1 to 1 transcripts per minute (6 to 60 transcripts per hour)
- Moderately transcribed genes: 1 to 10 transcripts per minute (60 to 600 transcripts per hour)
- Highly transcribed genes: greater than 10 transcripts per minute (over 600 transcripts per hour)


- **Rapidly degraded proteins:**  Degradation rate: >1 h^-1  Half-life: <1 hour
- **Stably expressed proteins:**   Degradation rate: between 0.1 and 1 h^-1   Half-life: anywhere from 1 to 10 hours
- **Very stable proteins:**   Degradation rate: <0.1 h^-1   Half-life: >10 hours

- **mRNA Degradation rate:** on the order of 0.5 to 3 min^-1 (30 to 180 h^-1)

 */

let prod_R3 = 3.0; // (3.33 minutes^-1 acc to Zhuo Chen) 1 (0.84 minutes^-1 acc to Shristopher A. Voigt) (0.84 acc to Jennifer S. Hallinan)
let prod_R1 = 3.0; //(3.33 minutes^-1 acc to Zhuo Chen) 1 (48 minutes^-1 acc to Shristopher A. Voigt)
let prod_L = 3.0; // (3.33 minutes^-1 acc to Zhuo Chen) 1 (6 minutes^-1 acc to Jennifer S. Hallinan)

let prod_MP1 = 2.0 ; // (1.416 transcript/minute acc to Zhuo Chen)1  (9 transcript/minute acc to Shristopher A. Voigt) (1.2 trancript/minute acc to Jennifer S. Hallinan)
let prod_MP3 = 2.0; // (1.667 transcript/minute acc to Zhuo Chen)1 (16.8 transcript/minute acc to Shristopher A. Voigt)
let prod_ML = 2.0; //   (2.083 transcript/minute acc to Zhuo Chen)1 (12 transcript/minute acc to Jennifer S. Hallinan)

let deg_R = 5e-3;  //  (3.33e-3 1/minute acc to Zhuo Chen) 1 (0.12 1/minute acc to Shristopher A. Voigt)
let deg_I = 5e-3; //  (3.33e-3 1/minute acc to Zhuo Chen) 1 (1.2 1/minute acc to Shristopher A. Voigt)
let deg_L = 0.015; // (0.01 1/minute acc to Zhuo Chen) 1

let deg_MP = 0.5; //  (0.138 1/minute acc to Zhuo Chen)1 

let for_com_RI =  9e-3; //  (5.33e-3 molecules^-1*min^-1 acc to Zhuo Chen)1 (9.06e-3 molecules^-1*min^-1 acc to Joseph A. Newman)
let for_com_RL = 2e-3; //  (5.33e-3 molecules^-1*min^-1 acc to Zhuo Chen)1 (1.61e-3 molecules^-1*min^-1 acc to Joseph A. Newman)

let n_R = 2; 
let n_A = 2; 
let K_R = 400; 
let K_A = 100; 

let init_I=0.005; // at the beginning, there is no SinI. AbrB is inhibiting the transcription of SinI and there's no Spo0A-P to activate it.
let init_L=21.8; // at the beginning, there is no SlrR. SinR is inhibiting the transcription of SlrR.
let init_m=0.05;

var initialValues = [];
for (var i = 0; i <= 20; i++) {
    initialValues.push([2380 , init_I, init_L, init_m, 4.0, 0.11, i*10]);
}
//  [R    , I      , L    , mI  , mR  , mL , A0 ];
var data1 = [];
var data2 = [];

for (var j = 0; j <  initialValues.length   ; j++) {

    var R = [initialValues[j][0]];
    var I = [initialValues[j][1]];
    var L = [initialValues[j][2]];
    var MP1 = [initialValues[j][3]];
    var MP3 = [initialValues[j][4]];
    var ML = [initialValues[j][5]];
    var P = [initialValues[j][6]];
    var c = [];

    for(i = 0; i < pointCount; i++) {

         A =P; 
        
        let activation_P1 = 1 / (1 + Math.pow(R[i] / K_R, n_R)) * Math.pow(A / K_A, n_A) / (1 + Math.pow(A / K_A, n_A));
       
        let activation_L = 1 / (1 + Math.pow(R[i] / K_R, n_R));
        // dMP3/dt equation

        let dMP3 = prod_MP3 - deg_MP*MP3[i];
        MP3.push(MP3[i] + dMP3/30);

        // dM1/dt equation
        let dMP1 = prod_MP1*activation_P1 - deg_MP*MP1[i];
        MP1.push(MP1[i] + dMP1/30);

        // dML/dt equation
        let dML = prod_ML*activation_L - deg_MP*ML[i];
        ML.push(ML[i] + dML/30);

        // dR/dt equation
         let dR = ((MP1[i])/10 + MP3[i])*prod_R3 - deg_R * R[i]- for_com_RI * R[i] * I[i];// - for_com_RL * R[i] * L[i]; /* + dis_com_L * ComL[i] + dis_com_I * ComI[i]; */
        R.push(R[i] + dR/30);

        
        // dI/dt equation
        let dI = (MP1[i])*prod_R1 - deg_I * I[i]- for_com_RI* R[i] * I[i];           
        // ; /* + dis_com_I * ComI[i]; */
        I.push(I[i] + dI/30);
        
        // dL/dt equation
        
        let dL = (ML[i])*prod_L - deg_L * L[i];
        L.push(L[i] + dL/30);

        c.push(i)
    }

    data1.push({
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
            cmax: pointCount
        },
        marker: {
            size: 1,
            color: c,
            colorscale: "Viridis",
            cmin: 0,
            cmax: pointCount
        }
    });

data2.push({
    type: 'scatter3d',
    mode: 'lines+markers',
    x: MP3, //SinR
    y: MP1, //SinI
    z: ML, //SlrR
    line: {
        width: 5,
        color: c,
        colorscale: "Viridis",
        cmin: 0,
        cmax: pointCount
    },
    marker: {
        size: 1,
        color: c,
        colorscale: "Viridis",
        cmin: 0,
        cmax: pointCount
    }
});

}
Plotly.newPlot('myDiv2', data2, {
    scene: {
        xaxis: {
            title: 'SinR mRNA [molecules]',
            titlefont: {
                color: 'red'
            }
        },
        yaxis: {
            title: 'SinI mRNA [molecules]',
            titlefont: {
                color: 'blue'
            }
        },
        zaxis: {
            title: 'SlrR mRNA [molecules]',
            titlefont: {
                color: 'rgb(255, 200, 0)'
            }
        }
    },
    showlegend: false
});

Plotly.newPlot('myDiv', data1, {
    scene: {
        xaxis: {
            title: 'SinR [molecules]',
            titlefont: {
                color: 'red'
            }
        },
        yaxis: {
            title: 'SinI [molecules]',
            titlefont: {
                color: 'blue'
            }
        },
        zaxis: {
            title: 'SlrR [molecules]',
            titlefont: {
                color: 'rgb(255, 200, 0)'
            }
        }
    },
    showlegend: false
});


