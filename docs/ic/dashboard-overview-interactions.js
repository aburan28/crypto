/* Presentation only: no fetch, solver dispatch, statistics, or run mutation. */
(() => {
  const library = document.getElementById('evidence-library');
  const ledger = document.getElementById('legacy-evidence');
  const query = document.getElementById('evidence-query');
  const list = document.getElementById('evidence-results');
  const count = document.getElementById('evidence-count');
  const reports = Array.from(ledger.querySelectorAll('h2')).map((heading, i) => {
    let target = heading.closest('[id]');
    if (!target || target === ledger) { heading.id = 'evidence-record-' + i; target = heading; }
    return { title: heading.textContent.trim(), id: target.id,
      text: (heading.closest('section') || heading.parentElement).textContent.toLowerCase() };
  });
  function search() {
    const words = query.value.toLowerCase().trim().split(/\s+/).filter(Boolean);
    const matches = reports.filter(report => words.every(word => report.text.includes(word)));
    list.replaceChildren();
    // Keep the first screen quiet until the reader asks for a report.
    if (!words.length) {
      count.textContent = reports.length + ' historical reports. Search by solver, curve or experiment.';
      return;
    }
    matches.slice(0, 10).forEach(report => {
      const li = document.createElement('li'); const a = document.createElement('a');
      a.href = '#' + report.id; a.textContent = report.title;
      a.addEventListener('click', () => { library.open = true; });
      li.append(a); list.append(li);
    });
    count.textContent = matches.length ? matches.length + ' matching reports' +
      (matches.length > 10 ? ' · showing 10; narrow your search' : '') :
      'No matching reports. Try a broader term, or open the full ledger.';
  }
  function revealHash() {
    let id;
    try { id = decodeURIComponent(location.hash.slice(1)); } catch { return; }
    const target = document.getElementById(id);
    if (!target) return;
    // Open every enclosing detail, including moved historical summaries.
    let parent = target.parentElement;
    while (parent) { if (parent.tagName === 'DETAILS') parent.open = true; parent = parent.parentElement; }
    requestAnimationFrame(() => target.scrollIntoView());
  }
  query.addEventListener('input', search);
  window.addEventListener('hashchange', revealHash);
  search(); revealHash();
})();
