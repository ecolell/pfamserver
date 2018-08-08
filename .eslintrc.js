module.exports = exports = {
    "extends": "standard",
    "env": {
      "browser": true,
      "amd": true,
      "jasmine": true
    },
    "globals": {
      "CONFIG": true,
      "$": true
    },
    "rules": {
      "comma-dangle": [
        "error", "always-multiline"
      ],
      "space-before-function-paren": [
        "error", "never"
      ],
      "semi": [
        "error", "always"
      ],
      "indent": [
        "error",
        4
      ],
      "quote": [
        "single"
      ],
      "eqeqeq": [
        "off"
      ],
      "camelcase": [
        "off"
      ],
      "no-mixed-operators": [
        "off"
      ],
      "no-unused-vars": [
        "warn",
        {
          "args": "all"
        }
      ],
      "one-var": [
        "off"
      ],
      "no-console": "warn",
      "handle-callback-err": [
        "off"
      ],
      "import/no-webpack-loader-syntax": [
        "off"
      ]
    }
}
