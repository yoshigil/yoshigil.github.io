/* ============================================================
   main.js — Small interactive behaviours
   JavaScript lets you respond to user actions (clicks,
   scrolls, etc.) and change the page dynamically.
   ============================================================ */

// --- 1. MOBILE NAV TOGGLE ---
// querySelector finds an element on the page by its CSS selector.
// addEventListener waits for an event (like a click) to happen.
const navToggle = document.querySelector('.nav__toggle');
const navLinks  = document.querySelector('.nav__links');

if (navToggle && navLinks) {
    navToggle.addEventListener('click', () => {
        // classList.toggle adds the class if missing, removes it if present
        navLinks.classList.toggle('open');
    });
}

// --- 2. HIGHLIGHT ACTIVE NAV LINK ---
// This reads the current page URL and adds an "active" class
// to whichever nav link matches, so it gets the underline styling.
const currentPath = window.location.pathname;
document.querySelectorAll('.nav__links a').forEach(link => {
    // link.getAttribute('href') gets the href value of the <a> tag
    const href = link.getAttribute('href');
    if (
        (href === 'index.html' && (currentPath === '/' || currentPath.endsWith('index.html'))) ||
        (href !== 'index.html' && currentPath.endsWith(href))
    ) {
        link.classList.add('active');
    }
});

// --- 3. FADE-IN ANIMATION ON SCROLL ---
// IntersectionObserver watches elements and fires a callback
// when they enter or leave the visible part of the screen.
const observer = new IntersectionObserver((entries) => {
    entries.forEach(entry => {
        if (entry.isIntersecting) {
            entry.target.classList.add('visible');
        }
    });
}, { threshold: 0.1 }); // fire when 10% of element is visible

// Add the observer to all cards and sections
document.querySelectorAll('.card, .pub-item, .highlight-box').forEach(el => {
    el.classList.add('fade-in');
    observer.observe(el);
});
