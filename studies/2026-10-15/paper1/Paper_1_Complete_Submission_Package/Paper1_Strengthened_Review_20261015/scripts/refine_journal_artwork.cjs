// Regenerate RGB PNGs from final SVGs after enforcing the publisher's 1 pt
// minimum line width. SVG text remains editable; Arial is supplied by macOS.
const fs = require('fs');
const path = require('path');
const runtime = process.env.P1_RUNTIME_MODULES || '/Users/zain/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules';
const sharp = require(path.join(runtime, 'sharp'));
const root = path.resolve(__dirname, '..');
const dir = path.join(root, 'publication_figures');
(async () => {
  const report = [];
  const manifest = JSON.parse(fs.readFileSync(path.join(root, 'results/publication_figure_manifest.json')));
  const current = new Set(manifest.map(x => x.filename + '.svg').concat(['GRAPHICAL_ABSTRACT.svg']));
  // Remove stale numbering copies from the final-artwork folder only. The
  // scientific source figures and pre-revision archive preserve older assets.
  for (const file of fs.readdirSync(dir).filter(x => /\.(svg|png)$/.test(x))) {
    if (!current.has(file.replace(/\.png$/, '.svg'))) fs.unlinkSync(path.join(dir, file));
  }
  for (const file of fs.readdirSync(dir).filter(x => x.endsWith('.svg'))) {
    const source = path.join(dir, file);
    let svg = fs.readFileSync(source, 'utf8');
    let changed = 0;
    const tag = svg.match(/<svg\b[^>]*>/)[0];
    const widthPt = Number(tag.match(/\bwidth="([0-9.]+)pt"/)[1]);
    const heightPt = Number(tag.match(/\bheight="([0-9.]+)pt"/)[1]);
    const intendedWidthPt = Math.min(6.9 * 72, 6.65 * 72 * widthPt / heightPt);
    const sourceScale = widthPt / intendedWidthPt;
    const minFont = Math.ceil(8 * sourceScale * 10) / 10;
    const minLine = Math.ceil(sourceScale * 100) / 100;
    svg = svg.replace(/stroke-width:\s*([0-9.]+)(?=[;"\s])/g, (match, number) => {
      const width = Number(number);
      if (width > 0 && width < minLine) { changed++; return 'stroke-width: ' + minLine; }
      return match;
    });
    // Matplotlib's default .8 pt stroke can be implicit on axes/tick paths.
    svg = svg.replace('stroke-linejoin: round; stroke-linecap: butt',
                      'stroke-linejoin: round; stroke-linecap: butt; stroke-width: ' + minLine);
    svg = svg.replace(/(font(?:-size)?:[^;]*?)([0-9.]+)px/g,
      (match, prefix, number) => prefix + Math.max(Number(number), minFont) + 'px');
    svg = svg.replace(/<text\b([^>]*)>([a-i]  [^<]*)<\/text>/g,
      (match, attrs, label) => '<text' + attrs.replace(/style="/, 'style="font-weight: bold; ') + '>' + label + '</text>');
    fs.writeFileSync(source, svg);
    const output = source.replace(/\.svg$/, '.png');
    // Explicit pixel dimensions avoid double density scaling in the bundled
    // SVG renderer. Preserve the original point-sized vector file.
    const rasterSvg = svg.replace(/<svg\b[^>]*>/, tag =>
      tag.replace(/\b(width|height)="([0-9.]+)pt"/g,
        (match, dimension, value) => dimension + '="' + Math.round(Number(value) * 600 / 72) + 'px"'));
    await sharp(Buffer.from(rasterSvg), {density: 72})
      .flatten({background: '#ffffff'}).toColourspace('srgb')
      .withMetadata({density: 600}).png().toFile(output);
    const meta = await sharp(output).metadata();
    report.push({filename: file, raised_strokes: changed, width: meta.width,
                 height: meta.height, density: meta.density, colourspace: meta.space,
                 intended_width_inches: intendedWidthPt / 72,
                 minimum_print_font_points: minFont / sourceScale,
                 minimum_print_stroke_points: minLine / sourceScale});
  }
  fs.writeFileSync(path.join(root, 'results/journal_artwork_validation.json'),
                   JSON.stringify(report, null, 2) + '\n');
  console.log('Finalized', report.length, 'RGB raster/vector figures at 600 dpi.');
})();
