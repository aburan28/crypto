# Cairn Security website

Static consultancy website for Cairn Security, covering post-quantum readiness,
migration strategy, implementation and integration, and cryptographic review.

## Preview locally

From the repository root:

```sh
python3 -m http.server 8000 --bind 127.0.0.1 --directory cairn-security
```

Open <http://localhost:8000>. No build step or package installation is required.
Relative asset paths also support serving the site from a subdirectory.

## Editing

- `index.html`: page copy, services, standards references, navigation, and favicon.
- `styles.css`: typography, responsive layouts, and charcoal/orange theme.
- `script.js`: mobile navigation and copyright year.
- `assets/cairn-hero.webp`: original AI-generated cairn artwork, optimized for the web.

The stylesheet requests DM Sans and Space Grotesk from Google Fonts; system fonts
are used if that service is unavailable. All other page assets are local.

The contact section deliberately says “Consultation contact details coming soon.”
Replace that text in `#contact-slot` with the approved business email or booking
link before accepting client enquiries. There is no enquiry form or backend.

## Hosted preview

The original private preview is hosted at
<https://cairn-security.aburan28.chatgpt.site> and requires owner access.
This directory preserves its design and content, with portable relative asset
paths. Commits here do not automatically update that hosted preview. No hosting
credentials or account-specific deployment manifest are included.

The website is independent of the educational Rust cryptography implementation
in this repository and does not execute any of that library's algorithms.

## Checks

```sh
node --check cairn-security/script.js
```

For content changes, check local asset references and section anchors. Before a
public launch, review the page at desktop and mobile widths and exercise the menu,
service disclosures, and keyboard navigation.
