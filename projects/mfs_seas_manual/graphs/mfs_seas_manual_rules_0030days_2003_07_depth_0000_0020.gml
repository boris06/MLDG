<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
<!-- Created by igraph -->
  <key id="name" for="graph" attr.name="name" attr.type="string"/>
  <key id="long_name" for="node" attr.name="long_name" attr.type="string"/>
  <key id="lat_wmc" for="node" attr.name="lat_wmc" attr.type="double"/>
  <key id="name" for="node" attr.name="name" attr.type="string"/>
  <key id="lon_wmc" for="node" attr.name="lon_wmc" attr.type="double"/>
  <key id="supp" for="edge" attr.name="supp" attr.type="double"/>
  <key id="conf" for="edge" attr.name="conf" attr.type="double"/>
  <graph id="G" edgedefault="directed">
    <data key="name">200307</data>
    <node id="n0">
      <data key="long_name">Alboran</data>
      <data key="lat_wmc">36.2056</data>
      <data key="name">ALB</data>
      <data key="lon_wmc">-2.6653</data>
    </node>
    <node id="n1">
      <data key="long_name">South-Western</data>
      <data key="lat_wmc">37.9649</data>
      <data key="name">SWE</data>
      <data key="lon_wmc">4.8654</data>
    </node>
    <node id="n2">
      <data key="long_name">North-Western</data>
      <data key="lat_wmc">41.4344</data>
      <data key="name">NWE</data>
      <data key="lon_wmc">5.044</data>
    </node>
    <node id="n3">
      <data key="long_name">Tyrrhenian</data>
      <data key="lat_wmc">39.9394</data>
      <data key="name">TYR</data>
      <data key="lon_wmc">12.0242</data>
    </node>
    <node id="n4">
      <data key="long_name">Adriatic</data>
      <data key="lat_wmc">42.9317</data>
      <data key="name">ADR</data>
      <data key="lon_wmc">15.8992</data>
    </node>
    <node id="n5">
      <data key="long_name">Ionian</data>
      <data key="lat_wmc">37.8056</data>
      <data key="name">ION</data>
      <data key="lon_wmc">18.9733</data>
    </node>
    <node id="n6">
      <data key="long_name">Central</data>
      <data key="lat_wmc">34.1235</data>
      <data key="name">CEN</data>
      <data key="lon_wmc">16.5527</data>
    </node>
    <node id="n7">
      <data key="long_name">Aegean</data>
      <data key="lat_wmc">38.0253</data>
      <data key="name">AEG</data>
      <data key="lon_wmc">24.9051</data>
    </node>
    <node id="n8">
      <data key="long_name">North-Levantine</data>
      <data key="lat_wmc">35.4952</data>
      <data key="name">NLE</data>
      <data key="lon_wmc">31.871</data>
    </node>
    <node id="n9">
      <data key="long_name">South-Levantine</data>
      <data key="lat_wmc">32.947</data>
      <data key="name">SLE</data>
      <data key="lon_wmc">28.6565</data>
    </node>
    <edge source="n4" target="n4">
      <data key="supp">0.00054</data>
      <data key="conf">0.9875</data>
    </edge>
    <edge source="n4" target="n5">
      <data key="supp">7e-006</data>
      <data key="conf">0.0125</data>
    </edge>
    <edge source="n7" target="n7">
      <data key="supp">0.000185</data>
      <data key="conf">0.771429</data>
    </edge>
    <edge source="n7" target="n6">
      <data key="supp">1.4e-005</data>
      <data key="conf">0.057143</data>
    </edge>
    <edge source="n7" target="n5">
      <data key="supp">2.7e-005</data>
      <data key="conf">0.114286</data>
    </edge>
    <edge source="n7" target="n9">
      <data key="supp">1.4e-005</data>
      <data key="conf">0.057143</data>
    </edge>
    <edge source="n0" target="n0">
      <data key="supp">3.4e-005</data>
      <data key="conf">0.384615</data>
    </edge>
    <edge source="n0" target="n2">
      <data key="supp">2.7e-005</data>
      <data key="conf">0.307692</data>
    </edge>
    <edge source="n0" target="n1">
      <data key="supp">2.7e-005</data>
      <data key="conf">0.307692</data>
    </edge>
    <edge source="n6" target="n6">
      <data key="supp">0.001731</data>
      <data key="conf">0.910072</data>
    </edge>
    <edge source="n6" target="n5">
      <data key="supp">7.5e-005</data>
      <data key="conf">0.039568</data>
    </edge>
    <edge source="n6" target="n9">
      <data key="supp">9.6e-005</data>
      <data key="conf">0.05036</data>
    </edge>
    <edge source="n5" target="n4">
      <data key="supp">2.1e-005</data>
      <data key="conf">0.042857</data>
    </edge>
    <edge source="n5" target="n6">
      <data key="supp">0.000137</data>
      <data key="conf">0.285714</data>
    </edge>
    <edge source="n5" target="n5">
      <data key="supp">0.000321</data>
      <data key="conf">0.671429</data>
    </edge>
    <edge source="n8" target="n7">
      <data key="supp">7e-006</data>
      <data key="conf">0.018182</data>
    </edge>
    <edge source="n8" target="n6">
      <data key="supp">2.1e-005</data>
      <data key="conf">0.054545</data>
    </edge>
    <edge source="n8" target="n5">
      <data key="supp">1.4e-005</data>
      <data key="conf">0.036364</data>
    </edge>
    <edge source="n8" target="n8">
      <data key="supp">0.000198</data>
      <data key="conf">0.527273</data>
    </edge>
    <edge source="n8" target="n9">
      <data key="supp">0.000137</data>
      <data key="conf">0.363636</data>
    </edge>
    <edge source="n2" target="n2">
      <data key="supp">0.000766</data>
      <data key="conf">0.817518</data>
    </edge>
    <edge source="n2" target="n1">
      <data key="supp">0.000116</data>
      <data key="conf">0.124088</data>
    </edge>
    <edge source="n2" target="n3">
      <data key="supp">5.5e-005</data>
      <data key="conf">0.058394</data>
    </edge>
    <edge source="n9" target="n6">
      <data key="supp">4.1e-005</data>
      <data key="conf">0.031088</data>
    </edge>
    <edge source="n9" target="n8">
      <data key="supp">0.000198</data>
      <data key="conf">0.150259</data>
    </edge>
    <edge source="n9" target="n9">
      <data key="supp">0.001081</data>
      <data key="conf">0.818653</data>
    </edge>
    <edge source="n1" target="n0">
      <data key="supp">7e-006</data>
      <data key="conf">0.011236</data>
    </edge>
    <edge source="n1" target="n6">
      <data key="supp">4.8e-005</data>
      <data key="conf">0.078652</data>
    </edge>
    <edge source="n1" target="n2">
      <data key="supp">0.000103</data>
      <data key="conf">0.168539</data>
    </edge>
    <edge source="n1" target="n1">
      <data key="supp">0.00039</data>
      <data key="conf">0.640449</data>
    </edge>
    <edge source="n1" target="n3">
      <data key="supp">6.2e-005</data>
      <data key="conf">0.101124</data>
    </edge>
    <edge source="n3" target="n6">
      <data key="supp">9.6e-005</data>
      <data key="conf">0.096552</data>
    </edge>
    <edge source="n3" target="n2">
      <data key="supp">7e-006</data>
      <data key="conf">0.006897</data>
    </edge>
    <edge source="n3" target="n1">
      <data key="supp">6.2e-005</data>
      <data key="conf">0.062069</data>
    </edge>
    <edge source="n3" target="n3">
      <data key="supp">0.000828</data>
      <data key="conf">0.834483</data>
    </edge>
  </graph>
</graphml>
