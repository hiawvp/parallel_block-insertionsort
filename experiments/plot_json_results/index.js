//@tscheck
const fs = require('fs');
const { ChartJSNodeCanvas } = require('chartjs-node-canvas');
const commandLineArgs = require('command-line-args')

function parseExtraArgs(extraArgs) {
  let parsedExtraArgs = "";
  if (extraArgs && extraArgs.length > 0) {
    parsedExtraArgs = "(";
    for (let index = 0; index < extraArgs.length; index++) {
      const element = extraArgs[index];
      parsedExtraArgs += `${element.name}=${element.value}`
      parsedExtraArgs += index < extraArgs.length - 1 ? ", " : "";
    }
    parsedExtraArgs += ")";
  }
  return parsedExtraArgs;
}

const optionDefinitions = [
  { name: 'input', alias: 'i', type: String, defaultOption: true },
  { name: 'output', alias: 'o', type: String },
  { name: 'plotTitle', alias: 't', type: String }
]
const options = commandLineArgs(optionDefinitions);

const inputFilename = options.input

if (inputFilename === undefined) {
  console.log("Missing required arg 'input'");
  process.exit(1);
}

if (!fs.existsSync(inputFilename)) {
  console.log(`input file ${inputFilename} not found`);
  process.exit(1);
}

let outputName;
if (options.output === undefined) {
  outputName = inputFilename.split('/').at(-1).split('.')[0];
} else {
  outputName = options.output;
}


// Load input JSON
// const inputFilename = '../cpu_ver/output_normal_1684277529.json'
const jsonData = fs.readFileSync(inputFilename);
const data = JSON.parse(jsonData);
const outputFormat = 'pdf';
console.log(inputFilename, outputName + "." + outputFormat);


let title;
if (options.plotTitle === undefined) {
  console.log(data.results[0].exec_times);
  let expInfo = data.experiment_info;
  let dataType = expInfo.input_type;
  let distribution = expInfo.input_distribution;
  let extraArgs = parseExtraArgs(expInfo.extra_args);
  let reps = data.results[0].exec_times[0].times.length;
  title = `Experiment ${dataType} ${distribution} ${extraArgs}, ${reps}  reps`;
} else {
  title = options.plotTitle;
}


const chartJSNodeCanvas = new ChartJSNodeCanvas({
  type: outputFormat,
  width: 800,
  height: 600,
  chartCallback: (ChartJS) => {
    // Global config example: https://www.chartjs.org/docs/latest/configuration/
    ChartJS.defaults.responsive = false;
    ChartJS.defaults.maintainAspectRatio = false;
    ChartJS.defaults.elements.line.borderWidth = 2;
  }
});

// Define the chart options...
const pointStyles = ['circle'
  , 'cross', 'crossRot', 'dash', 'line', 'rect', 'rectRounded', 'rectRot', 'star', 'triangle']

const colors = ['#CB4335', '#1F618D', '#F1C40F', '#27AE60', '#884EA0', '#D35400'];

// Create a new chart for each element in the input JSON array
const datasets = data.results.map((element, idx) => {
  // Extract data from the current element
  const name = element.name;
  const execTimes = element.exec_times;


  const borderDash = idx % 2 == 0 ? undefined : [5, 5];
  const pointStyle = pointStyles[Math.floor(Math.random() * pointStyles.length)];

  // Create a new dataset
  const dataset = {
    label: name + parseExtraArgs(element.extra_args),
    data: [],
    fill: false,
    showLine: true,
    borderDash,
    pointStyle,
    pointRadius: 5,
    borderColor: colors[idx % colors.length],
  };

  // Add data points to the dataset
  execTimes.forEach((time) => {
    const n = time.n;
    const times = time.times;
    if (times.length > 0) {
      const avgTime = times.reduce((acc, val) => acc + val) / times.length;
      dataset.data.push({ x: n, y: avgTime });
    }
  });

  // Add the dataset to the chart options
  // chartOptions.data.datasets.push(dataset);
  return dataset;
});

const lines = [];
// const lines = datasets.map((element) => {
//   element.type = 'line';
//   return element;
// })
// console.log(lines);


/** @type {import("chart.js").ChartOptions} */
const chartOptions = {
  type: 'scatter',
  data: {
    datasets: [...datasets, ...lines],
  },
  options: {
    plugins: {
      title: {
        display: true,
        text: title
      },
      legend: {
        labels: {
          usePointStyle: true,
        },
      }
    },
    scales: {
      x: {
        min: 0,
        // bounds: 'data',
        // type: 'logarithmic',
        // offset: true,
        // type: 'timeseries',
        title: {
          display: true,
          text: 'Input size (n)'
        },
        ticks: {
          // align: 'end',
          suggestedMin: 1000,
          maxTicksLimit: 5
        },
        gridLines: {
          display: false
        }
      },
      y: {
        type: 'logarithmic',
        display: true,
        ticks: {
          padding: 10,
        },
        title: {
          display: true,
          text: 'Execution time (s)'
        },
      }
    },
    title: {
      display: true,
      text: 'Execution time vs. input size'
    }
  }
};


// Generate the chart and write it to a file
(async () => {
  const canvas = chartJSNodeCanvas.renderToBufferSync(chartOptions);
  fs.writeFileSync(`${outputName}.pdf`, canvas);
})();

