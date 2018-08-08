'use strict';

// imports
const webpack = require('webpack');
const path = require('path');
const ExtractTextPlugin = require('extract-text-webpack-plugin');
const ManifestRevisionPlugin = require('manifest-revision-webpack-plugin');
const CleanWebpackPlugin = require('clean-webpack-plugin');


// paths
const rootAssetPath = path.resolve(__dirname, './web');

// env variables
const RUN_ENV = {
    DEV: 'dev',
    PRODUCTION: 'production',
    WEBPACK_DEV_SERVER: 'webpack-dev-server'
};

const UglifyJSPlugin = require('uglifyjs-webpack-plugin');
const MergeFilesPlugin = require('merge-files-webpack-plugin');
const PrepackWebpackPlugin = require('prepack-webpack-plugin').default;
const SitemapPlugin = require('sitemap-webpack-plugin').default;


const dotenv = require('dotenv');
const appPath = rootAssetPath + '/src/client/apps/';
const nodeEnv = process.env.NODE_ENV === undefined ? RUN_ENV.DEV : process.env.NODE_ENV;

function parseChunkValue(name, value) {
    let ext = value.split('.').slice(-1)[0];
    return name + '.' + ext;
}

let myCoolFormatter = function (data, parsedAssets) {
    let output = {};
    output.assets = parsedAssets;

    for (let chunkName in data.assetsByChunkName) {
        let chunkValue = data.assetsByChunkName[chunkName];
        if (typeof chunkValue === 'string') {
            output.assets[parseChunkValue(chunkName, chunkValue)] = chunkValue;
        } else if (Array.isArray(chunkValue)) {
            chunkValue.forEach(function (val) {
                output.assets[parseChunkValue(chunkName, val)] = val;
            });
        }
    }
    output.publicPath = data.publicPath;
    return output;
};


const paths = [
  '#!/job_request/',
  '#!/job_result/example1',
  '#!/job_result/example2',
  '#!/job_result/example3',
  '#!/help',
];


module.exports = {
    devServer: {
        inline: false,
        headers: {
            'Access-Control-Allow-Origin': '*'
        },
        compress: true
    },
    entry: {
        main: [
            appPath + 'main.js'
        ]
    },
    output: {
        path: __dirname + '/pfamserver/static/dist',
        publicPath: 'http://localhost:2992/assets/',
        filename: nodeEnv === RUN_ENV.DEV || nodeEnv === RUN_ENV.WEBPACK_DEV_SERVER
            ? '[name].[chunkhash:8].js'
            : '[name].[chunkhash:8].min.js',
        chunkFilename: '[id].[chunkhash:8].js'
    },
    resolve: {
        extensions: ['*', '.js', '.html', '.css'],
    },
    resolveLoader: {
        moduleExtensions: ["-loader"]
    },
    module: {
        loaders: [
            {
                enforce: 'pre',
                test: /\.js$/,
                exclude: [/node_modules/, /shared/, /vendor/, /test/],
                loader: "eslint-loader", // todo: replace for 'standard-loader', once the all fixes are completed.
                options: {
                    emitError: true,
                    emitWarning: true,
                    failOnError: true,
                }
            },
            {
                test: /\.(jpe?g|png|gif)$/i,
                loader: "file-loader",
                options: {
                    name: 'img/[name].[hash:8].[ext]',
                }
            },
            {
                test: /\.js$/i,
                exclude: [/node_modules/, /test/],
                loader: 'babel-loader',
                query: {
                    presets: ['env']
                }
            },
            {
                test: /\.(css|less)$/,
                use: ExtractTextPlugin.extract({
                    fallback: "style-loader",
                    use: [{
                        loader: "css-loader",
                        options: {
                            minimize: true
                        }
                    },
                    {
                        loader: "less-loader",
                        options: {
                            minimize: true
                        }
                    }]
                })
            },
            {
                test: /\.scss$/,
                use: ExtractTextPlugin.extract({
                    fallback: 'style-loader',
                    use: [{
                        loader: "css-loader",
                        options: {
                            minimize: true
                        }
                    },
                        'sass-loader'
                    ]
                }),
                exclude: [/node_modules/, /test/]
            },
            {
                test: /\.html$/,
                loader: "html",
                exclude: [/node_modules/, /test/]
            },
            {
                test: /\.(woff2|woff|ttf|eot|svg)(\?v=[a-z0-9]\.[a-z0-9]\.[a-z0-9])?$/,
                loader: "url-loader"
            },
            {
                test: /\.py$/,
                use: 'raw-loader'
            }
        ]
    },
    plugins: [
        new webpack.LoaderOptionsPlugin({
            test: /\.js$/,
            exclude: (nodeEnv === RUN_ENV.PRODUCTION ? [/test/] : undefined),
            options: {
              eslint: {
                configFile: '.eslintrc.js', // this is my helper for resolving paths
                cache: false,
              }
            },
        }),
        //new webpack.optimize.CommonsChunkPlugin({ name: 'main', filename: appPath + 'main.js' }),
        new ExtractTextPlugin({
            filename: '[name].[hash:6].css',
            allChunks: true
        }),
        new ManifestRevisionPlugin(path.join('pfamserver/static/', 'manifest.json'), {
            rootAssetPath: rootAssetPath,
            ignorePaths: ['/src', '/test'],
            format: myCoolFormatter
        }),
        new webpack.ProvidePlugin({
            jQuery: 'jquery',
            $: 'jquery',
            jquery: 'jquery',
            ko: 'knockout'
        }),
        new MergeFilesPlugin({
            filename: 'css/style.css',
            test: /min\.css/,
            deleteSourceFiles: true
        }),
        //new PrepackWebpackPlugin({})
        new SitemapPlugin('https://pfamserver.leloir.org.ar', paths, {
            changeFreq: 'daily',
            priority: '0.4'
        }),
  ]
};


//Configure /dist folder for non nn server environments
if (nodeEnv !== RUN_ENV.WEBPACK_DEV_SERVER) {
    module.exports.plugins.push(new webpack.optimize.UglifyJsPlugin({
        compress: {
            warnings: false
        },
        output: {
            comments: false
        }
    }));
    module.exports.output.publicPath = '/static/dist/';
    module.exports.plugins.push(new CleanWebpackPlugin(['pfamserver/static/dist/'], {
        verbose: true,
        dry: false
    }));
}

//Configure Prod build
if (nodeEnv === RUN_ENV.PRODUCTION) {
} else {
    //Configure non-Prod build
    module.exports.plugins.push(new webpack.DefinePlugin({
        'CONFIG': JSON.stringify(
            {
            })
    }));
    module.exports.devtool = 'cheap-source-map';
}

module.exports.node = {
    fs: 'empty',
    child_process: 'empty',
    net: 'empty'
};

