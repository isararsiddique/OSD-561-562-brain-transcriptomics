// MathJax 3.2.2 display equations; install mathjax-full when not bundled.
const fs=require('fs'); const path=require('path');
const modules=process.env.P1_MATH_MODULES||'/tmp/p1_equation_renderer/node_modules';
const runtime=process.env.P1_RUNTIME_MODULES||'/Users/zain/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const {mathjax}=require(path.join(modules,'mathjax-full/js/mathjax.js'));
const {TeX}=require(path.join(modules,'mathjax-full/js/input/tex.js'));
const {SVG}=require(path.join(modules,'mathjax-full/js/output/svg.js'));
const {liteAdaptor}=require(path.join(modules,'mathjax-full/js/adaptors/liteAdaptor.js'));
const {RegisterHTMLHandler}=require(path.join(modules,'mathjax-full/js/handlers/html.js'));
const {AllPackages}=require(path.join(modules,'mathjax-full/js/input/tex/AllPackages.js'));
const sharp=require(path.join(runtime,'sharp'));
const adaptor=liteAdaptor();RegisterHTMLHandler(adaptor);
const doc=mathjax.document('',{InputJax:new TeX({packages:AllPackages}),OutputJax:new SVG({fontCache:'local'})});
const equations=[
String.raw`A_p=I-vv^T,\quad s=w^TY_0,\quad\gamma=\frac{Y_0s^T}{ss^T},\quad A_s=I-\gamma w^T.\tag{1}`,
String.raw`A_s^2=A_s,\quad w^TA_s=0,\quad\|A_s\|_2=\|\gamma\|_2=\sqrt{1+\|h\|_2^2}.\tag{2}`,
String.raw`c=\min\left\{1,\frac{\sqrt{L^2-1}}{\|h\|_2}\right\},\quad\gamma_L=w+ch,\quad A_L=I-\gamma_L w^T.\tag{3}`,
String.raw`\lambda_{\max}(A_L\Sigma A_L^T)\le L^2\lambda_{\max}(\Sigma),\quad\|A_Lb-b\|_2\le L|w^Tb|.\tag{4}`,
String.raw`\beta=b+\delta u\quad\Longrightarrow\quad A_L\beta-b=-\gamma_L(w^Tb)+\delta\{u-\gamma_L(w^Tu)\}.\tag{6}`,
String.raw`\mathbb E\|\widehat\gamma\|^2=1+\frac{G-1}{r-2},\qquad r>2.\tag{7}`,
String.raw`\operatorname{Var}(\|\widehat\gamma\|^2)=\frac{2(G-1)(G+r-3)}{(r-2)^2(r-4)},\qquad r>4.\tag{8}`,
String.raw`\operatorname{Energy\ removal}(S)=100\left[1-\frac{\sum_{g\in S}(\beta_g^{*})^2}{\sum_{g\in S}\beta_g^2}\right].\tag{5}`,
String.raw`C^*=\underset{\|C\|_2\le\sqrt{L^2-1}}{\arg\min}\left\{\frac{\|Y-SC\|_F^2}{2n}+\frac{\alpha}{2}\|DC\|_F^2\right\}.\tag{9}`,
String.raw`A=\begin{pmatrix}0&0\\-C^T&I_q\end{pmatrix},\quad A^2=A,\quad E_sA=0,\quad\|A\|_2=\sqrt{1+\|C\|_2^2}\le L.\tag{10}`,
String.raw`\|Ab-b\|_2\le L\|E_sb\|_2.\tag{11}`,
String.raw`\operatorname{gap}(T)=\langle\nabla f(T),T\rangle+\sqrt{L^2-1}\|\nabla f(T)\|_*\ge f(T)-f(T^*).\tag{12}`];
(async()=>{const root=path.resolve(__dirname,'..');const dir=path.join(root,'figures/equations');fs.mkdirSync(dir,{recursive:true});let info=[];
const numberingPath=path.join(root,'results/strengthening_equation_number_map.json');
if(fs.existsSync(numberingPath)){const numbering=JSON.parse(fs.readFileSync(numberingPath,'utf8'));for(let i=0;i<equations.length;i++)equations[i]=equations[i].replace(/\\tag\{(\d+)\}/g,(_,n)=>`\\tag{${numbering[n]}}`);}
for(let i=0;i<equations.length;i++){
 let h=adaptor.outerHTML(doc.convert(equations[i].replace(/\\tag\{([^}]+)\}/g, '\\qquad\\text{($1)}'),{display:true}));let svg=h.slice(h.indexOf('<svg'),h.indexOf('</svg>')+6);
 if(!svg.includes('xmlns='))svg=svg.replace('<svg ','<svg xmlns="http://www.w3.org/2000/svg" ');
 svg=svg.replace(/(width|height)="([0-9.]+)(ex|em)"/g,(_,a,n)=>`${a}="${Math.round(Number(n)*8.5)}px"`).replace(/currentColor/g,'#000000');
 if(!svg.match(/width="([0-9.]+)px"/)){ console.log(svg.slice(0,600)); throw new Error('Responsive SVG needs explicit sizing'); }
 const width=Number(svg.match(/width="([0-9.]+)px"/)[1]);const name='equation_'+(i+1);
 fs.writeFileSync(path.join(dir,name+'.svg'),svg);await sharp(Buffer.from(svg),{density:600}).flatten({background:'#ffffff'}).png().toFile(path.join(dir,name+'.png'));
 info.push({filename:name,width_inches:Math.min(6.85,width/96),latex:equations[i]});
}fs.writeFileSync(path.join(dir,'manifest.json'),JSON.stringify(info,null,2)+'\n');console.log('Rendered',info.length,'display equations.');})();
