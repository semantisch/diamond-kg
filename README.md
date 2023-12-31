## Badges
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)

<img src="diamond-kg.png" width=250 height=250 />

# Diamond-KG

A toolbox for constructing Knowledge Graphs from natural language sentences from drug indications
with focus on their medical context and therapeutic intent.

## Features

- Extract Chat-GPT-generated triples
- Extract medical context (free and defined context types)
- Identify links with NIH NCATS Translator SRI Name Resolution API

## Requirements

- Python 3.6+
- openai 0.27.8+
- rdflib 6.2.0+
- Requests 2.31.0+
- rich 13.4.2+
- urllib3 1.26.13+

## Installation

You can clone this repository to your machine and install the required packages from requirements.txt.

Make sure the following requirements are fulfilled:
- [You have an OpenAI key with access to the gpt-4 model.](https://openai.com/gpt-4)

## Usage

```bash
  diamond-kg.py 
    --input <path to the input JSON file> 
    --prompt ["triples" | "freeContext" | "definedContext"] 
    --output <path to the output JSON file> 
    --apiKey <OpenAI API key with gpt-4 model access> 
    --organization <OpenAI API organization>
    [--chunkSize <Number of paragraphs per prompt>]
    [--maxTokens <Max tokens>]
```

## Input

For the structure of input, refer to neurodkg_sentences.json.

## Result datasets

The example result datasets can be found in the results directory of this repository.

## Contributing

Contributions are welcome. Please open an issue or submit a pull request if you would like to contribute.

## License

DIAMOND-KG is licensed under the [MIT](https://opensource.org/licenses/MIT) license.
