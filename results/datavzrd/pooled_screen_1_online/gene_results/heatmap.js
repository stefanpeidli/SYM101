const heatmap_plot = {
    "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
    "data": {
        "values": []
    },
    "hconcat": [
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "[object]",
                    "align":"left"
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'[object]'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "text": {
                        "field": "",
                        "type": "[object]",
                        "title": "",
                        "legend": null
                    }
                }
            }
            
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Gscreend_FDR'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Gscreend_FDR",
                        "type": "quantitative",
                        "title": "Gscreend_FDR",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Gscreend_FDR"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Gscreend_pvalue'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Gscreend_pvalue",
                        "type": "quantitative",
                        "title": "Gscreend_pvalue",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Gscreend_pvalue"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Gscreend_LFC'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Gscreend_LFC",
                        "type": "quantitative",
                        "title": "Gscreend_LFC",
                        "scale": {
                            "domain": ["-5","5"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Gscreend_LFC"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Gscreend_cell_line_specificity'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Gscreend_cell_line_specificity",
                        "type": "quantitative",
                        "title": "Gscreend_cell_line_specificity",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "oranges",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Gscreend_cell_line_specificity"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_filtered_score'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_filtered_score",
                        "type": "quantitative",
                        "title": "MAGeCK_filtered_score",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "oranges",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_filtered_score"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_filtered_p-value'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_filtered_p-value",
                        "type": "quantitative",
                        "title": "MAGeCK_filtered_p-value",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_filtered_p-value"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_filtered_fdr'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_filtered_fdr",
                        "type": "quantitative",
                        "title": "MAGeCK_filtered_fdr",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_filtered_fdr"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_filtered_lfc'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_filtered_lfc",
                        "type": "quantitative",
                        "title": "MAGeCK_filtered_lfc",
                        "scale": {
                            "domain": ["-5","5"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_filtered_lfc"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_score'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_score",
                        "type": "quantitative",
                        "title": "MAGeCK_score",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "oranges",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_score"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_p-value'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_p-value",
                        "type": "quantitative",
                        "title": "MAGeCK_p-value",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_p-value"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_fdr'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_fdr",
                        "type": "quantitative",
                        "title": "MAGeCK_fdr",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_fdr"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'MAGeCK_lfc'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "MAGeCK_lfc",
                        "type": "quantitative",
                        "title": "MAGeCK_lfc",
                        "scale": {
                            "domain": ["-5","5"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "MAGeCK_lfc"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Bagel2_filtered_Bayesfactor'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Bagel2_filtered_Bayesfactor",
                        "type": "quantitative",
                        "title": "Bagel2_filtered_Bayesfactor",
                        "scale": {
                            "domain": ["500","2000"],
                            "type": "linear",
                            "scheme": "oranges",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Bagel2_filtered_Bayesfactor"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Bagel2_filtered_FDR'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Bagel2_filtered_FDR",
                        "type": "quantitative",
                        "title": "Bagel2_filtered_FDR",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Bagel2_filtered_FDR"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Bagel2_Bayesfactor'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Bagel2_Bayesfactor",
                        "type": "quantitative",
                        "title": "Bagel2_Bayesfactor",
                        "scale": {
                            "domain": ["500","2000"],
                            "type": "linear",
                            "scheme": "oranges",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Bagel2_Bayesfactor"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]},
        
        {
            
            "layer": [
            
            {
                "mark": {
                    "type": "rect",
                    "align":"left",
                    "width": 20
                },
                "encoding": {
                    "x": {
                        "field": "dummy",
                        "axis": {
                            "labelExpr": "'Bagel2_FDR'",
                            "labelAngle": -90,
                            
                            "title": null,
                            "ticks": false,
                            "orient": "top",
                            "domain": false
                        }
                    },
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    },
                    "color": {
                        "field": "Bagel2_FDR",
                        "type": "quantitative",
                        "title": "Bagel2_FDR",
                        "scale": {
                            "domain": ["0","1"],
                            "type": "linear",
                            "scheme": "blues",
                            
                            "clamp": true,
                        },
                        "legend": null
                    }
                }
            }
            ,
            {
                "mark": {
                    "type": "text",
                    "align":"left",
                    "xOffset": -6
                },
                "encoding": {
                    "text": {"field": "Bagel2_FDR"},
                    "y": {
                        "field": "index",
                        "type": "ordinal",
                        "axis": null,
                    }
                }
            }
        
        ]}
        
    ],
    "config": {
        "style": {"cell": {"stroke": "transparent"}, "guide-label": {"fontWeight": "bold"}},
        "concat": {"spacing": 0},
        "text": {"limit": 135}
    }
}